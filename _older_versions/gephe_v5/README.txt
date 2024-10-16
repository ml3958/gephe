---------
Author: Menghan Liu (ml464@columbia.edu)
Saeed Tavazoie lab, Department of Biology, Columbia University
---------


---------
Usage:
---------

For alignment

+---------------+
| Determine if  |
| folder exists |
+---------------+
         |
+-----------+----------+
|                      |
+------+------+      +--------+--------+
| Folder exists? |      | Folder is empty? |
+------+------+      +--------+--------+
|                      |
+---+---+              +---+---+
| Yes    |              | No    |
|        |              |      |
| Prompt |              |  X   |
+---+----+              +---+---+
|                      |
+------+-------+      +-------+-------+
| User responds |      |  X function is  |
| with input    |      |   performed     |
+------+-------+      +----------------+
|
+------+-------+
| End program  |
+--------------+


---------
File structure:
---------
- dir_diamond:
      - input.fasta
      - input_database.dmnd
      - 649633097.diamond.out
      - 649633097.diamond.out.pickle

- dir_pp:
      - top_*%_*.csv_alignment.mci
      - top_*%_*.csv_alignment.mci_I*.mcl
      - top_*%_*.csv_alignment.mci_I*.mcl_annot_full.txt
      - top_*%_*.csv_alignment.mci_I*.mcl_annot_summary.txt
      - top_*%_*.csv_alignment.mci_I*.mcl.pickle
      - top_*%_*.csv_alignment.mci_I*.mcl_pp_ortho.csv
      - top_*%_*.csv.dic
      - top_*%_*.csv_pp.csv
      - test

- dir_pp_pog:
