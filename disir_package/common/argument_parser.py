import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description="DiSiR is a computational strategy for "
                                                 "identifying ligand-receptor interactions at cell type "
                                                 "level from scRNA-seq data.")

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-sp", "--scRNA-path", help="Path to scRNA-Seq data", type=str,
                          default=None, required=True)
    required.add_argument("-ctp", "--cell-type-path", help="Path to cell type labels", type=str,
                          default=None, required=True)
    required.add_argument("-gp", "--gene-path", help="Path to gene names", type=str,
                          default=None, required=True)
    required.add_argument("-sip", "--subunit-interactions-path", help="Path to ligand-receptor interactions at subunit level", type=str,
                          default=None, required=True)
    required.add_argument("-odp", "--output-directory-path", help="Path to output directory where the output results will be saved", type=str,
                          default=None, required=True)

    # I/O options
    parser.add_argument("-iv", "--iteration-value", help="Number of iterations for permutating data in statistical test", type=int,
                        default=1000)
    parser.add_argument("-tf", "--threshold-fraction", help="Threshold on fraction of cells expressed each ligand or receptor per cell type", type=float,
                        default=0)
    parser.add_argument("-te", "--threshold-expression", help="Threshold on scaled (max-normalized) average expression of each ligand or receptor within a cell type", type=float, 
                        default=0)
    parser.add_argument("-tp", "--threshold-pvalue", help="Threshold on p-value for filtering non-significant LR interactions", type=float, 
                        default=0.05)
    parser.add_argument("-sm", "--sparse-matrix", help="Use sparse scRNA matrix (for big input scRNA data) in the Harwell-Boeing or MatrixMarket (.mtx) format", action="store_true")

    arguments = parser.parse_args()

    return arguments.__dict__
