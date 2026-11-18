"""Verify both polynomial rings. With --all, regenerate rows and test finite rings."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

def main():
    if not __debug__: raise RuntimeError("Run with assertions enabled")
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--all",action="store_true")
    args=ap.parse_args();here=Path(__file__).resolve().parent
    for line in (here/"SHA256SUMS").read_text().splitlines():
        checksum,name=line.split("  ",1)
        if hashlib.sha256((here/name).read_bytes()).hexdigest()!=checksum:
            raise ValueError("checksum mismatch: "+name)
    def call(script,*argv):
        return subprocess.check_output([sys.executable,str(here/script),*map(str,argv)],cwd=here)
    call("check_bezout.py","--verify",here/"bezout_weights.json")
    result=[json.loads(s) for s in call("verify_rows.py",*(here/("rows%d.json"%n) for n in (6,7,8))).splitlines()]
    if result!=json.loads((here/"row_checks.json").read_bytes()): raise ValueError("row checks")
    print("PASS: both polynomial rings at n=6,7,8 and all 16 Bezout identities.",flush=True)
    with tempfile.TemporaryDirectory() as td:
        target=Path(td)/"direct.json"
        call("direct_quotient.py",here/"bezout_weights.json",*(here/("rows%d.json"%n) for n in (6,7,8)),"--output",target)
        if json.loads(target.read_bytes())!=json.loads((here/"direct_checks.json").read_bytes()):
            raise ValueError("direct quotient")
    print("PASS: division-free quotient in both rings at orders 6, 7 and 8.",flush=True)
    if not args.all: return
    with tempfile.TemporaryDirectory() as td:
        target=Path(td)/"colours.json"
        call("check_colour_trace.py","--max-strands","6","--write",target)
        if json.loads(target.read_bytes())!=json.loads((here/"colour_trace.json").read_bytes()):
            raise ValueError("shaded projector normalization")
    print("PASS: independent two-colour trace and denominator checks through six strands.",flush=True)
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td)
        for n in (6,7,8):
            new=tmp/("rows%d.json"%n)
            call("construct_rows.py",n,new)
            if new.read_bytes()!=(here/new.name).read_bytes(): raise ValueError(new.name)
        call("verify_certificate.py",*(here/name for name in
             ("certificate4.json","colour_certificate4.json","certificate5.json","colour_certificate5.json")))
        for orders,name in [("1,2,3,4,5,6,7","ramified_fibres.json"),("8","ramified_fibres_order8.json")]:
            new=tmp/name
            call("check_ramified_fibres.py","--orders",orders,"--write",new)
            if json.loads(new.read_bytes())["checks"]!=json.loads((here/name).read_bytes())["checks"]:
                raise ValueError(name)
    print("PASS: generators reproduced, small explicit matrices and 16 finite-ring checks.")

if __name__=="__main__":
    main()
