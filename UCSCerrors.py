from SSerrors import load_errors,write_errors,print_errors,Error

class DataFileError(Error):

    error_type = "DataFileError"
    def __init__(self, code, message):
        self.code = code
        self.message = message

class SequenceDataError(Error):
    error_type = "SequenceDataError"
    def __init__(self, code, message):
        self.code = code
        self.message = message

class NCBIQueryError(Error):
    error_type = "NCBIQueryError"
    def __init__(self, code, message):
        self.code = code
        self.message = message
