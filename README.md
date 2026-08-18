# HIV-Mpox-Research

- Log 1 (July 19, 2026): Remade every Phase I code, ready to run.
- Log 2 (July 19, 2026): Ran each Phase I code successfully, took 2 hours to complete.
- Log 3 (August 14, 2026): Started remaking Phase II code.
- Log 4 (August 15, 2026): Completed a mock prototype of Phase II, revisions needed to be applied in actual work but should work as intended.
- Log 5 (August 16, 2026): Completed Phase II steps A and B, ready to run, made Github repo to store the codes and its results.
- Log 6 (August 16, 2026): Completed Phase II step C, changed MolPro for DeepSol and ran AggreScan3D locally instead.
- Log 7 (August 17, 2026): Completed Phase II step D, marking the finished product for Phase II but reliability remains a problem as ramachandran flavored plot yielded 94.85 which is less than the 95 threshold expected. The group decided to round it up and proceed instead, opting to re-do phase I and II if phase III yields bad results.
- Log 8 (August 18, 2026): Re did phase II step D entirely with the integration of phenix and was able to produce a usable pdb file  that hit the targeted >95% ramachandran flavored percentage and fixed structure.


ADDITIONAL INFO:
- All runs on Python and was executed locally.
- For Phase I, work as is yan dapat unless path finding is hard wired to my Mac.
- Phase II required Biopython and pdb2pqr which I installed in a virtual environment alongside FoldX (nakuha ko sa website nila ung binary) and freesa which was installed locally on my device through terminal.
- Phase II Common file is a requirement for all Phase II code. Make sure to create a folder "_common" and paste the phase2_common file inside the containing folder of all codes.
- Phase II flows follows:
    Execute Phase II A -> Execute Phase II B -> Import Sequences to AlphaFold 3 -> Import Files to Phase II B -> Execute Phase II B -> Execute Phase II C -> Import Files to DeepSol, AggreScan3D, and CamSol -> Import Files to Phase II C -> Execute Phase II C -> Execute Phase II D -> Import files to MolProb -> Import Files to Phase II D -> Execute Phase II D. 
