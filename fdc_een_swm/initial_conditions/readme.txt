Run:

python3 polvani_conditions.py Rossby Froude Nx outfilestem --csv

    The (optional) CSV flag will give you 4 CSV files: outfilestem_u.txt, outfilestem_v.txt, outfilestem_h.txt, outfilestem_data.txt which contain the u, v, and h in array order (transpose to real-space) and a data file with some things like Rossby, Froude, Nx, etc. Omitting this flag gives you an ugly bundle of stuff which plays nicely with my code instead.

    The domain size, periodicity, scaling, etc. are the same as before. Please let me know if you have any questions or concerns. For a domain size of 256, I find it takes between 10-30 iterations to converge. Every now and again it seems to get stuck and just oscillate around, so if this happens just kill it and try again.
    
Case B
python3 polvani_conditions.py 0.05 0.05 256 outfilestem --csv
