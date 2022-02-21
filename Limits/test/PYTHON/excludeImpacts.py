import json
from copy import deepcopy
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="input JSON file")
parser.add_argument("-o", "--output", dest="output", type=str, default="", help="output JSON file (if empty: append '_include' to input file name)")
parser.add_argument("-m", "--matches", dest="matches", type=str, default=[], nargs='*', help="exclude nuisances that match args")
parser.add_argument("-x", "--exact-matches", dest="exact_matches", type=str, default=[], nargs='*', help="exclude nuisances that exactly match args")
args = parser.parse_args()

if len(args.output)==0:
    args.output = args.input.replace(".json","_include.json")

with open(args.input,'r') as infile:
    input = json.load(infile)
output = deepcopy(input)

filtered_params = []
for param in input["params"]:
    if any(m in param["name"] for m in args.matches) or any(x==param["name"] for x in args.exact_matches):
        continue
    else:
        filtered_params.append(param)

output["params"] = filtered_params

# dump in the format used by combineTool/Impacts
with open(args.output,'w') as outfile:
    json.dump(output, outfile, sort_keys=True, indent=2, separators=(',', ': '))
