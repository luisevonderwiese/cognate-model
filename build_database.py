from lingdata import database

config_path = "prototype_lingdata_config.json"
database.read_config(config_path)
#database.download()
database.compile()
