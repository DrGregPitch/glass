import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from fastapi.testclient import TestClient
import app as appmod
client = TestClient(appmod.app)

def test_elements_endpoint():
    j = client.get("/api/elements").json()
    els = j["elements"]; assert len(els) == 118 and "PubChem" in j["source"]
    fe = els[25]; assert fe["symbol"] == "Fe" and fe["row"] == 4 and fe["col"] == 8 and abs(fe["mass"] - 55.84) < 0.01
    assert {(e["row"], e["col"]) for e in els}.__len__() == 118          # unique layout cells
    assert els[56]["row"] == 9 and els[88]["row"] == 10                    # La / Ac series rows
