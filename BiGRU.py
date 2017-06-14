import os
os.environ['KERAS_BACKEND'] = 'theano'
from keras.models import model_from_yaml
import numpy as np




class GRUNetwork():
    
    def __init__(self, model_name = 'BiGRU.yaml', weights_name = "BiGRU.h5", mem_depth = 20):
        yaml_file = open(model_name, 'r')
        loaded_model_yaml = yaml_file.read()
        yaml_file.close()
        self.model = model_from_yaml(loaded_model_yaml)
        self.model.load_weights(weights_name)
        self.mem_depth = 20
        self.alphabet_structures = ['C', 'H', 'E', 'T', "*"]
        self.alphabet_proteins = ['A','R','N','D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        
        self.structures_indices = dict((c, i) for i, c in enumerate(self.alphabet_structures))
        self.indices_structures = dict((i, c) for i, c in enumerate(self.alphabet_structures))

        self.proteins_indices = dict((c, i) for i, c in enumerate(self.alphabet_proteins))
        self.indices_proteins = dict((i, c) for i, c in enumerate(self.alphabet_proteins))
        
    def secondary_structure(self, seqq):
        
        seqq_with_borders = np.append(np.append(["*"]*int(self.mem_depth/2 - 1), list(seqq)), ["*"]*int(self.mem_depth/2))
        
        #Get time series
        protein_blocks = []
        for i in range(0, len(seqq_with_borders) - self.mem_depth + 1):
            protein_blocks.append(seqq_with_borders[i: i + self.mem_depth])


        #Vectorisation
        X = np.zeros((len(protein_blocks), self.mem_depth, len(self.alphabet_proteins)), dtype=np.bool)
        for i, block in enumerate(protein_blocks):
            for t, protein in enumerate(block):
                X[i, t, self.proteins_indices[protein]] = 1
                
        
        
        predictions = (self.model.predict(X))
        predictions = [self.alphabet_structures[np.argmax(prediction)] for prediction in predictions]
        
        return ''.join(predictions)