"""Verify the supplement and regenerate the data in a temporary directory."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

def main():
    if not __debug__: raise RuntimeError("Run with assertions enabled")
    here=Path(__file__).resolve().parent
    for line in (here/"SHA256SUMS").read_text().splitlines():
        checksum,name=line.split("  ",1)
        if hashlib.sha256((here/name).read_bytes()).hexdigest()!=checksum:
            raise ValueError("checksum mismatch: "+name)
    def call(script,*args):
        return subprocess.check_output([sys.executable,str(here/script),*map(str,args)],cwd=here)
    with tempfile.TemporaryDirectory() as td:
        temp=Path(td);joint=temp/"joint.json";longest=temp/"longest.json"
        call("joint_counts.py",here/"b005315.txt",joint)
        call("longest_runs.py",joint,longest)
        for fresh,frozen in [(joint,"joint_counts.json"),(longest,"longest_runs.json")]:
            if json.loads(fresh.read_bytes())!=json.loads((here/frozen).read_bytes()):
                raise ValueError(frozen)
        current=json.loads(call("verify_surgery.py"))
        expected=json.loads((here/"surgery_checks.json").read_bytes())
        for key in ("parameters","orders","script_sha256"):
            if current[key]!=expected[key]: raise ValueError("surgery "+key)
    with tempfile.TemporaryDirectory() as td:
        target=Path(td)/"positivity.json"
        call("check_positivity.py",here/"b005315.txt",here/"joint_counts.json","--output",target)
        if json.loads(target.read_bytes())!=json.loads((here/"positivity_checks.json").read_bytes()):
            raise ValueError("shifted positivity")
        for script,parameters,frozen in [
            ("check_idle_moves.py",[here/"b005315.txt","--output",target],"idle_moves_checks.json"),
            ("check_return_words.py",["--write",target],"return_words.json"),
            ("weighted_check.py",["--joint",here/"joint_counts.json","--write",target],"weighted_checks.json"),
        ]:
            call(script,*parameters)
            if json.loads(target.read_bytes())!=json.loads((here/frozen).read_bytes()):
                raise ValueError(frozen)
    print("PASS: sharp shift, all 741 differences, independent return-word counts and weighted transfers.")
    print("PASS: both histograms through 38 and all exhaustive surgery checks through 8.")

if __name__=="__main__":
    main()
