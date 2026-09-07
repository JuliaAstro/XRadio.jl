#!/usr/bin/env python3
"""Fetch the reference MSv4 processing set the UVW-convention test reads.

The store is a public S3 bucket of many small Zarr chunks rather than a single
file, so it is walked key by key rather than fetched with wget.
"""

import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BUCKET = "https://viper-test-data.s3.amazonaws.com/"
STORE = "Antennae_North.cal.lsrk.split.v8.ps.zarr"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def keys(prefix):
    token = None
    while True:
        query = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            query["continuation-token"] = token
        url = BUCKET + "?" + urllib.parse.urlencode(query)
        root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
        for entry in root.findall(NS + "Contents"):
            yield entry.find(NS + "Key").text
        token = root.findtext(NS + "NextContinuationToken")
        if root.findtext(NS + "IsTruncated") != "true" or not token:
            return


def fetch(key, into):
    path = os.path.join(into, key)
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # The bucket resets connections under a burst, so back off and retry.
    for attempt in range(5):
        try:
            with urllib.request.urlopen(
                BUCKET + urllib.parse.quote(key), timeout=60
            ) as response, open(path, "wb") as out:
                out.write(response.read())
            return
        except OSError:
            time.sleep(1 + attempt)
    raise SystemExit(f"could not fetch {key}")


if __name__ == "__main__":
    into = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    for key in keys(STORE):
        fetch(key, into)
