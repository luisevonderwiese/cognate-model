from lingdata import database

config_path = "setup_comparison_lingdata_config.json"
database.read_config(config_path)
database.update_native()
database.generate_data()
