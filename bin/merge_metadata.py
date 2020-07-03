#!/usr/bin/env python3

import json
import sys


main_json = sys.argv[1]
metadata_json = sys.argv[2]
output_json = sys.argv[3]



if len(sys.argv) != 4:
    print("usage: {} <main_json> <metadata_json> <output_json>".format(
        sys.argv[0]))
    exit(1)


# Read metadata file
with open(metadata_json) as m:
    metadata = json.load(m)

# Read main report file
with open(main_json) as r:
    report = json.load(r)

# merge main report content with metadata
merge_report = metadata.copy()
merge_report.update(report)

# write final report to a file
with open(output_json, "w") as f:
    json.dump(merge_report, f, indent=4)
