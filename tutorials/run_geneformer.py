from geneformer import TranscriptomeTokenizer

# if previously wrote the 'ensembl_id', the following might not be necessary.
import scanpy as sc 
tem = sc.read_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/NSCLC_subsetted_raw.h5ad")
tem.var['ensembl_id'] = tem.var['features']
tem.__dict__['_raw'].__dict__['_var'] = tem.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'ensembl_id'})
tem.write_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/NSCLC_subsetted_raw.h5ad")


tk = TranscriptomeTokenizer({"cell_type": "cell_type"}, nproc=2)


tk.tokenize_data("/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/", 
                 "/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/", 
                 "NSCLC_subsetted", 
                 file_format="h5ad")


from geneformer import EmbExtractor
# initiate EmbExtractor
embex = EmbExtractor(model_type="CellClassifier",
                     num_classes=3,
                     filter_data={"cell_type":["Cardiomyocyte1","Cardiomyocyte2","Cardiomyocyte3"]},
                     max_ncells=None,
                     emb_layer=0,
                     emb_label=["disease","cell_type"],
                     labels_to_plot=["disease"],
                     forward_batch_size=200,
                     nproc=16)

# extracts embedding from input data
# example dataset: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/example_input_files/cell_classification/disease_classification/human_dcm_hcm_nf.dataset
embs = embex.extract_embs("../fine_tuned_models/geneformer-6L-30M_CellClassifier_cardiomyopathies_220224",
                          "path/to/input_data/",
                          "path/to/output_directory/",
                          "output_prefix")

