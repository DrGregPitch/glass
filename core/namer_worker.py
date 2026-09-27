"""
Sidecar process: SMILES -> candidate IUPAC names.  Runs inside namer_env (Python 3.11) because
the generator needs PyTorch/transformers.  Protocol: one SMILES per stdin line; one JSON line
per result on stdout: {"smiles":..., "candidates":[{"style":"SYST","name":...}, ...]}.
Generator: knowledgator/SMILES2IUPAC-canonical-base (MT5 seq2seq, Apache-2.0), asked for the
systematic, base and traditional styles so the OPSIN gate in core/namer.py can pick a verified one.
"""
import sys, json, os
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
try:
    from chemicalconverters import NamesConverter
    conv = NamesConverter(model_name="knowledgator/SMILES2IUPAC-canonical-base")
except Exception as e:  # pragma: no cover
    print(json.dumps({"error": f"generator unavailable: {e}"}), flush=True)
    sys.exit(1)
print(json.dumps({"ready": True, "generator": "knowledgator/SMILES2IUPAC-canonical-base"}), flush=True)
STYLES = ["SYST", "BASE", "TRAD"]
for line in sys.stdin:
    smi = line.strip()
    if not smi:
        continue
    try:
        cands, seen = [], set()
        for st in STYLES:                       # several beams per style: the OPSIN gate picks the first verified one
            outs = conv.smiles_to_iupac(f"<{st}>{smi}", num_beams=5, num_return_sequences=3)
            if isinstance(outs, str):
                outs = [outs]
            for n in outs:
                n = (n or "").strip()
                if n and n not in seen:
                    seen.add(n); cands.append({"style": st, "name": n})
        print(json.dumps({"smiles": smi, "candidates": cands}), flush=True)
    except Exception as e:
        print(json.dumps({"smiles": smi, "error": str(e)}), flush=True)
