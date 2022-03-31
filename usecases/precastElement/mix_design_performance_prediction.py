"""
predict properties of a mix"""

from argparse import ArgumentParser

def solve_and_write_output(
        water_cement_ratio: float,
    outputfile: str,
):
    """
    outputfile : str
        FilePath to the output file into which the CO2 emission is written.
    """
    # this should rather be a function of the cement content, but just for testing purposes
    CO2 = water_cement_ratio
    with open(outputfile, "w") as handle:
            handle.write("CO2:{}\n".format(CO2))


if __name__ == "__main__":
    PARSER = ArgumentParser(description="performance prediction of mix")
    PARSER.add_argument(
        "-wz",
        "--water_cement_ratio",
        required=True,
        help="water cement ratio",
    )
    PARSER.add_argument(
        "-o",
        "--outputfile",
        required=True,
        help="file name for the output to be written",
    )
    ARGS = vars(PARSER.parse_args())

    solve_and_write_output(ARGS["water_cement_ratio"], ARGS["outputfile"]
                           )