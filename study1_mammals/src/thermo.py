"""ThermoMPNN wrapper: per-site 20-value energy vectors for a backbone with an arbitrary threaded sequence.

ThermoMPNN predicts ddG(wt->mut) at site i as f_i[mut] - f_i[wt], where f_i (21 values) comes from frozen ProteinMPNN
embeddings of the structure and its sequence. Returning f_i for every site gives all ddG values in the context of
the threaded sequence with one forward pass. Positive ddG = destabilizing (ThermoMPNN convention).
"""
import sys

import numpy as np
import torch
from omegaconf import OmegaConf

from common import TOOLS

TM_DIR = TOOLS / 'ThermoMPNN-main'
sys.path.insert(0, str(TM_DIR))
from protein_mpnn_utils import alt_parse_PDB, tied_featurize  # noqa: E402
from transfer_model import TransferModel  # noqa: E402

TM_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
AA20 = 'ACDEFGHIKLMNPQRSTVWY'


def load_model(device=None):
    device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
    cfg = OmegaConf.create({'platform': {'thermompnn_dir': str(TM_DIR)},
                            'model': {'hidden_dims': [64, 32], 'subtract_mut': True, 'num_final_layers': 2,
                                      'freeze_weights': True, 'load_pretrained': True, 'lightattn': True}})
    model = TransferModel(cfg)
    state = torch.load(TM_DIR / 'models' / 'thermoMPNN_default.pt', map_location='cpu', weights_only=True)['state_dict']
    missing, unexpected = model.load_state_dict({k[len('model.'):]: v for k, v in state.items()}, strict=True)
    return model.eval().to(device)


def parse(pdb_path, chain='A'):
    return alt_parse_PDB(str(pdb_path), chain)[0]


@torch.no_grad()
def site_energies(model, pdb_dict, sequence=None):
    """Return array (L, 20) of f_i[a] in AA20 order; ddG(a->b at i) = E[i,b] - E[i,a]."""
    d = dict(pdb_dict)
    if sequence is not None:
        seq = sequence.replace('U', 'C')
        assert len(seq) == len(d['seq_chain_A'])
        d['seq_chain_A'] = seq
        d['seq'] = seq
    device = next(model.parameters()).device
    feats = tied_featurize([d], device, None, None, None, None, None, None, ca_only=False)
    X, S, mask, _, chain_M, chain_enc, *_rest = feats
    residue_idx = feats[12]
    hidden, embed, _ = model.prot_mpnn(X, S, mask, chain_M, residue_idx, chain_enc, None)
    h = torch.cat(hidden[:model.num_final_layers], -1)[0]          # (L, 256)
    x = torch.cat([h, embed[0]], -1).unsqueeze(-1)                  # (L, 384, 1): each site is its own length-1 input
    x = model.light_attention(x, mask)                             # (L, 384)
    out = model.ddg_out(model.both_out(x).unsqueeze(-1)).squeeze(-1)  # (L, 21)
    idx = [TM_ALPHABET.index(a) for a in AA20]
    return out[:, idx].float().cpu().numpy()


@torch.no_grad()
def reference_ddg(model, pdb_dict, positions):
    """Original per-mutation code path of ThermoMPNN, for validating site_energies."""
    from datasets import Mutation
    seq = pdb_dict['seq_chain_A']
    muts = [Mutation(position=i, wildtype=seq[i], mutation=b, ddG=None, pdb='x') for i in positions for b in AA20]
    pred, _ = model([pdb_dict], muts)
    return np.array([p['ddG'].item() for p in pred]).reshape(len(positions), 20)
