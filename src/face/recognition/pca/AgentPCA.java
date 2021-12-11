package face.recognition.pca;

/*-----------------------------------------------------------------------------
 File         : AgentPCA
  Author      : Engin DURMAZ
  Date        : 24.04.2019
  Description : Class which manages the training and test operations
               For the training phase 
                 - Calculate mean vector
                 - Extract features
                 - Apply projection
               For test phase
                 - Calculate distance
                 - find nearest neighbor
-----------------------------------------------------------------------------*/
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AgentPCA 
{
	private final int NEIGHBOR_K = 1;
	
	// trainingSet: label and image matrix pair which are trained
	private ArrayList<LabelMatrixPair> trainingSet = new ArrayList<>();
	// projectedSet: projected training set
	private ArrayList<LabelMatrixPair> projectedSet = null;
	// vector of mean values
	private Matrix meanMatrix = null;
	// vector of feature values
	private Matrix featureMatrix = null;
	// dimension of extraction
	private int dimensionOfExtraction = 3;
	
	/************** PRIVATE FUNCTIONS **********************/

	/*
	 * calculate means and fill mean vector
	 */
	private void calculateMean() 
	{
		int rowCount = trainingSet.get(0).getMatrix().getRowCount();
		int setCount = trainingSet.size();
		
		meanMatrix = new Matrix(rowCount, 1);

		/*
		 * visit all rows (column count is 1 in vector), calculate means
		 * 		1) sum of items in rows of all sets
		 * 		2) divide total to set count
		 * 		3) save results into mean vector
		 */
		for(int row = 0; row < rowCount; row++)
		{
			double rowTotal = 0;
			for (int trSet = 0; trSet < setCount; trSet++) 
			{
				rowTotal +=trainingSet.get(trSet).getMatrix().get(row, 0);
			}
			meanMatrix.set(row, 0, rowTotal/setCount);
		}
	}
	
	/*
	 * Extract features to a matrix
	 */
	private Matrix applyFeatureExtraction() 
	{
		int i, j;

		int row = trainingSet.get(0).getMatrix().getRowCount();
		int column = trainingSet.size();
		
		Matrix X = new Matrix(row, column);

		// firstly find differences from mean
		for (i = 0; i < column; i++) 
		{
			X.copyFrom(0, row - 1, i, i, 
					trainingSet.get(i).getMatrix().subtract(this.meanMatrix));
		}

		/*
		 *  get eigenvalues and eigenvectors
		 *  1. Transpose and 
		 *  2. Multiply by itself
		 */
		Matrix XT = X.transpose();
		Matrix XTX = XT.multiply(X);
		
		EigenOperations feature = new EigenOperations(XTX);
		double[] d = feature.getd();

		if( d.length < dimensionOfExtraction)
		{
			return null;// TODO: throw Exception: Not enough feature
		}
		
		int[] indexes = this.getIndexesOfKEigenvalues(d, dimensionOfExtraction);

		Matrix eigenVectors = X.multiply(feature.getV());
		Matrix selectedEigenVectors = getMatrix(eigenVectors, 0,
				eigenVectors.getRowCount() - 1, indexes);

		// normalize the eigenvectors
		row = selectedEigenVectors.getRowCount();
		column = selectedEigenVectors.getColumnCount();
		for (i = 0; i < column; i++) 
		{
			double temp = 0;
			for (j = 0; j < row; j++)
			{
				temp += Math.pow(selectedEigenVectors.get(j, i), 2);
			}
			
			temp = Math.sqrt(temp);

			for (j = 0; j < row; j++) 
			{
				selectedEigenVectors.set(j, i, selectedEigenVectors.get(j, i) / temp);
			}
		}

		return selectedEigenVectors;

	}

	private int[] getIndexesOfKEigenvalues(double[] eigenValues, int k) 
	{
		IndexValue [] ivPairs = new IndexValue[eigenValues.length];
		
		for (int i = 0; i < eigenValues.length; i++)
		{
			ivPairs[i] = new IndexValue(i, eigenValues[i]);
		}

		Arrays.sort(ivPairs);

		int[] result = new int[k];
		for (int i = 0; i < k; i++)
		{
			result[i] = ivPairs[i].index;
		}
		return result;
	}

	private Matrix getMatrix (Matrix srcMatrix, int rowStart, int rowEnd, int[] indices) 
	{
		Matrix target = new Matrix(rowEnd - rowStart + 1, indices.length);

		for (int i = rowStart; i <= rowEnd; i++) 
		{
			for (int j = 0; j < indices.length; j++) 
			{
				target.set(i - rowStart, j, srcMatrix.get(i, indices[j]));
			}
		}

		return target;
	}

	private int findMaxNeighbor(List<LabelMatrixPair> neighbors) 
	{
		int maxIndex = 0;
		for(int i = 1; i < NEIGHBOR_K; i++)
		{
			if(neighbors.get(i).getDistance() > neighbors.get(maxIndex).getDistance())
			{
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	/****************************** PUBLIC FUNCTIONS ****************************/
	
	public AgentPCA(int dimension)
	{
		dimensionOfExtraction = dimension;
	}
	
	public void addTrainingData(LabelMatrixPair pair)
	{
		trainingSet.add(pair);
	}
	/*
	 * Function applies PCA on training data finds projected set
         - Calculate means by calling calculateMean
         - cretae matrix of features by calling applyFeatureExtraction
         - Compute projection
	 */
	public void train()
	{
		calculateMean();
		featureMatrix = applyFeatureExtraction();
		
		projectedSet = new ArrayList<LabelMatrixPair>();
		for (int i = 0; i < trainingSet.size(); i++) 
		{
			/*
			 * Feature Extraction to training set
			 * W^T * (A(i) - mu(i))
			 */
			
			// 1. Subtract mean from current matrix
			Matrix diff = trainingSet.get(i).getMatrix().subtract(meanMatrix);
			
			// 2. Transpose Feature Matrix and multiply with diff matrix found
			Matrix extractedFeatures = featureMatrix.transpose().multiply(diff);
			LabelMatrixPair extractionPair = new LabelMatrixPair(extractedFeatures,
					trainingSet.get(i).getLabel());
			
			this.projectedSet.add(extractionPair);
		}
		
	}
	/*
	  findNearestNeighbor finds the most suitable matrix for
	  the given image matrix
         - Calculate distance
         - find nearest neighbor
	 */
	public LabelMatrixPair findNearestNeighbor(Matrix matrix) 
	{
		/*
		 * Feature Extraction to training set
		 * W^T * (A(i) - mu(i))
		 */
		
		// 1. Subtract mean from current matrix
		Matrix diff = matrix.subtract(meanMatrix);
		
		// 2. Transpose Feature Matrix and multiply with diff matrix found
		Matrix extractedFeatures = featureMatrix.transpose().multiply(diff);
		
		List<LabelMatrixPair> neighbors = new ArrayList<LabelMatrixPair>();
		/*
		 * Firstly assume first indices as nearest neighbors
		 */
		for(int i = 0; i < NEIGHBOR_K; i++)
		{
			LabelMatrixPair matrixPair = projectedSet.get(i); 
			double distance = projectedSet.get(i).getMatrix().getEuclideanDistance(extractedFeatures);
			matrixPair.setDistance(distance);
			neighbors.add(matrixPair);
		}
		/*
		 * Then, visit all items and compare visited item with found neighbors
		 */
		for(int i = NEIGHBOR_K; i < projectedSet.size(); i++)
		{
			LabelMatrixPair matrixPair = projectedSet.get(i); 
			double distance = projectedSet.get(i).getMatrix().getEuclideanDistance(extractedFeatures);
			matrixPair.setDistance(distance);
			
			int indexOfMax = findMaxNeighbor(neighbors);
			if(matrixPair.getDistance() < neighbors.get(indexOfMax).getDistance())
			{
				neighbors.remove(indexOfMax);
				neighbors.add(matrixPair);
			}
		}
        String labelFound = neighbors.get(0).getLabel();
        /*
         * Finally find maximum ratio of visited count per distance
         */
        for(int i = 0; i < NEIGHBOR_K; i++)
        {
        	
        }
        return neighbors.get(0);
    }
}
