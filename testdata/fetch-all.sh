#!/bin/sh

wget -r --reject "index.html*" -np -nH --cut-dirs=3 https://dl-dsoc.nrao.edu/ptiede/183766358/269cf6f450965db446809f78c0bbc9a3/BT164/BT164A///VLBA_BT164A_bt164aQband_BIN0_SRC0_0_260108T230453.idifits
wget -r --reject "index.html*" -np -nH --cut-dirs=3 https://dl-dsoc.nrao.edu/ptiede/183766358/269cf6f450965db446809f78c0bbc9a3/BT164/BT164A///VLBA_BT164A_bt164aWbandrecor_BIN0_SRC0_0_260113T200047.idifits

# The reference MSv4 store the UVW-convention test reads.
python3 "$(dirname "$0")/fetch-msv4.py" "$(dirname "$0")"
