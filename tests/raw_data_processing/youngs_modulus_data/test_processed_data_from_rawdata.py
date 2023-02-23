from lebedigital.raw_data_processing.youngs_modulus_data.emodul_generate_processed_data \
        import processed_data_from_rawdata
from pathlib import Path

def test_emodul_generate_processed_data():
        """
        Testing the processed data generation on a single experiment
        """

        # setup paths and directories
        data_dir = 'test_data'
        data_path = Path(__file__).parent / data_dir
        output_file = data_path / 'processed_data.csv'

        processed_data_from_rawdata(locationOfRawData=data_path, locationOfProcessedData=output_file)

        assert output_file.is_file

