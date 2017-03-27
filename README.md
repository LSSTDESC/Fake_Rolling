# Fake_Rolling
Small piece of code to fake a rolling cadence from a default Opsim output.
Merge fields of observations.

To run this code:

python Run_Fake_Rolling.py --fieldid 309 310 311 --dbFile /sps/lsst/data/dev/pgris/sims_operations/DB_Files/minion_1016_sqlite.db --merge_factor 0.5 --new_DB_name Rolling_minion_1016_309_310_311_80.db

where:
--fieldid : list of fields to be merged (9 at max)
--dbFile : Opsim file to be processed
--merge factor: fraction of fields to be merged

If only one field is mentioned, the merging is the following:
- season_1 = season_1
- season_2 = season2 + merge_factor*season_3 + merge_factor*season_4
- season_3 = (1.-merge_factor) * season_3
- season_4 = (1.-merge_factor) * season_4
and the same cycle for seasons 5 to 7 and 8 to 10

If there are multiple fields : illustration with three fields A, B and C:

- season_1_A = season_1_A
- season_1_B = season_1_B
- season_1_C = season_1_C

- season_2_A = season_2_A + merge_factor * season_2_B + merge_factor * season_2_C
- season_2_B = (1.-merge_factor) * season_2_B
- season_2_C = (1.-merge_factor) * season_2_C

- season_3_A = (1.-merge_factor) * season_3_A 
- season_3_B = season_3_B + merge_factor * season_3_A + merge_factor * season_3_C
- season_3_C = (1.-merge_factor) * season_3_C

- season_4_A = (1.-merge_factor) * season_4_A 
- season_4_B = (1.-merge_factor) * season_4_B
- season_4_C = season_4_C + merge_factor * season_4_A + merge_factor * season_4_B

and this cycle is made as long as it is allowed.