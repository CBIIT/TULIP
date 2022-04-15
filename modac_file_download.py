import utils
import os

model_url = "https://modac.cancer.gov/api/v2/dataObject/NCI_DOE_Archive/JDACS4C/JDACS4C_Pilot_1/Tumor_type_classifier_models/"

cnn_17_pc = "cnn_17_pc_weights.h5" 
cnn_17 = "cnn_17_weights.h5"
cnn_32_pc = "cnn_32_pc_weights.h5"
cnn_32 = "cnn_32_weights.h5"

models = [cnn_17_pc, cnn_17, cnn_32_pc, cnn_32]

for model in models:
    data_loc = utils.fetch_file(model_url + model, unpack=False, md5_hash=None, subdir="models")
    print('Data downloaded and stored at: ' + data_loc)
    data_path = os.path.dirname(data_loc)
    print(data_path)
