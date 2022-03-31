"""
solution precast element
"""

from argparse import ArgumentParser

def solve_and_write_output(
    CO2 : float,
    outputfile: str,
):
    """
    outputfile : str
        FilePath to the output file into which the CO2 emission is written.
    """
    volume = 10
    KPI_C02 = volume * CO2
    with open(outputfile, "w") as handle:
            handle.write("KPI_CO2:{}\n".format(CO2))


if __name__ == "__main__":
    PARSER = ArgumentParser(description="run structural simulation of precast element")
    PARSER.add_argument(
        "-o",
        "--outputfile",
        required=True,
        help="file name for the output to be written",
    )
    PARSER.add_argument(
        "-CO2",
        "--CO2",
        required=True,
        help="CO2 emission per m³ of concrete",
    )
    ARGS = vars(PARSER.parse_args())

    solve_and_write_output(ARGS["CO2"], ARGS["outputfile"]
                           )