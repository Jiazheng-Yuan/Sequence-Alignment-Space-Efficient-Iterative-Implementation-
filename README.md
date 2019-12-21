# Sequence-Alignment-Space-Efficient-Iterative-Implementation-
This is the project for CS466 Bio informatics. The project is the implementation for Hirschberg algorithm. Both iterative and recursive version
are implemented. The function entrance for recursive version is Hirschberg_Entrance in Alignment_Recursive.py. The function entrance for iterative version is 
Hirschberg in Alignment_Iterative.py. The default scoring function delta_fitting is provided at the top for both files. It is customizible.
The NeedlemanWunsch function is provided in both file, which is the space-non-efficient version. The evaluation function in both files calculate 
the alignment score, given the alignment, which can be use to verify that the alignment by NeedlemanWunsch and Hirschberg and Hirschberg_Entrance
has the same score, given the same input and delta.

The main function in both files are written for testing purpose. They can be called from commandline. Each takes two argument. which are the name of test files.
Test files should be stored in dataset folder, with the format seq1_size(size).txt and seq2_size(size).txt. (size) should be replaced with the actual length of the sequence in the file.
and this length should used as the commandline argument. size of seq1_size(size).txt will be used as the first commandline argument,size of seq2_size(size).txt will be used as the second commandline argument.
There are pre-generated testdata in dataset folder already. For example, python3 Alignment_Iterative 400 400 will align the sequences in seq1_size400.txt and seq2_size400.txt.


