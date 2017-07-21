
from Bio import SeqIO
#import sys
import os
import fire


class Annotation:

    def prediamond(self, root, output):
        """ Concatenates assembled contigs into one file so you can 

        Parameters:
        -----------
            root : str
                path to root directory where contigs are kept
            output : str
                path to the output file

        Examples:
        ---------
            >>> annon = Annotation()
            >>> annon.prediamond("/path/to/root", "/path/to/output/file")
        """
        with open(output, "w") as combined:
            kos = os.listdir(root)
            for ko in kos:
                file = "{}/{}/454AllContigs.fna".format(root, ko)
                with open(file, "r") as fa:
                    for entry in SeqIO.parse(fa, format="fasta"):
                        entry.id = "{}:{}".format(ko, entry.id)
                        SeqIO.write(entry, combined, format="fasta")

        bashCommand = """
        qsub -N diamond_contigs -V -cwd -b y -q all.q -pe orte 24 \\
            docker run --rm \\
            -v /scratch/uesu/:/scratch \\
            -v {}:/w/query.fna \\
            etheleon/diamond:0.1 diamond blastx --log -v -c 1 -b10.0 \\
            -q /w/query.fna -p 24 -t /scratch -d /scratch/db/nrfull.dmnd -o /scratch/output.m8 -f 6
        """.format(output)
        print(bashCommand)
if __name__ == '__main__':
    fire.Fire(Annotation)
