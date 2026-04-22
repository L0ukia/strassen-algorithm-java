public class Strassen {
   
	public static int strassOperations = 0;

	public static double[][] addMatrix(double[][] Amatrix,double[][] Βmatrix){
	// Πρόσθεση
		int n=Amatrix.length;
		double [][] C= new double[n][n];
		for(int i=0; i<n; i++ )
			for(int j=0; j<n; j++)
				C[i][j]=Amatrix[i][j]+Βmatrix[i][j];
		
		return C;
	}
	
	public static double[][]subMatrix(double[][] Amatrix,double[][] Bmatrix){
	// Αφαίρεση
		int n=Amatrix.length;
		double [][] C= new double[n][n];
		for(int i=0; i<n; i++ )
			for(int j=0; j<n; j++)
				C[i][j]=Amatrix[i][j]-Bmatrix[i][j];
		
		return C;
	}
	
	public static double[][]multMatrix(double[][] Amatrix,double[][] Bmatrix){
	//Πολλαπλασιασμός
		int n=Amatrix.length;
		double [][] C= new double[n][n];
		for(int i=0; i<n; i++ )
			for(int j=0; j<n; j++)
				C[i][j]=Amatrix[i][j]*Bmatrix[i][j];
		
		return C;
	}

	public static double[][] matrixMultiply(double[][] Amatrix, double[][] Bmatrix){
		
	        if(Amatrix.length != Bmatrix.length){
            throw new RuntimeException("Matrices not equal dimensions");
	        }
	        double[][] Cmatrix = new double[Amatrix.length][Bmatrix.length];
		int convOperations=0;
		for(int i = 0; i < Amatrix.length; i++){
	            for(int j = 0; j < Amatrix.length; j++){
	                for(int k = 0; k < Amatrix.length; k++){
	                    Cmatrix[i][j] += Amatrix[i][k]*Bmatrix[k][j];
	                    convOperations += 2; // μια πράξη πολλαπλασιασμού και μια πρόσθεση
	                }
	            }
		}
		return Cmatrix;
	    }

	public static double[][] extractSubMatrix(double[][] matrix, int row, int col, int size) {
	  //Εξαγωγή υποπίνακα
		double[][] submatrix = new double[size][size];
	    for (int i = 0; i < size; i++) {
	        for (int j = 0; j < size; j++) {
	            submatrix[i][j] = matrix[row + i][col + j];
	        }
	    }
	    return submatrix;
	}

	public static double[][] strassen(double[][] firstMat, double[][] secondMat) {
		
        if (firstMat.length != secondMat.length) {
            throw new RuntimeException("Matrices not equal dimensions");
        }
        int n=firstMat.length;
        if (n == 1) {
            firstMat[0][0] = firstMat[0][0] * secondMat[0][0];
            strassOperations++;
            return firstMat;
        }
      
        double[][] productMat = new double[n][n];

       // Μέρος 1
        double[][] a11 = extractSubMatrix(firstMat, 0, 0, n / 2);
        double[][] a12 = extractSubMatrix(firstMat, 0, n / 2, n / 2);
        double[][] a21 = extractSubMatrix(firstMat, n / 2, 0, n / 2);
        double[][] a22 = extractSubMatrix(firstMat, n / 2, n / 2, n / 2);

        double[][] b11 = extractSubMatrix(secondMat, 0, 0, n / 2);
        double[][] b12 = extractSubMatrix(secondMat, 0, n / 2, n / 2);
        double[][] b21 = extractSubMatrix(secondMat, n / 2, 0, n / 2);
        double[][] b22 = extractSubMatrix(secondMat, n / 2, n / 2, n / 2);

        int operationCost = a11.length * a11.length;

       // Μέρος 2
        double[][] s1 = subMatrix(b12, b22);
        strassOperations += operationCost;
        double[][] s2 = addMatrix(a11, a12);
        strassOperations += operationCost;
        double[][] s3 = addMatrix(a21, a22);
        strassOperations += operationCost;
        double[][] s4 = subMatrix(b21, b11);
        strassOperations += operationCost;
        double[][] s5 = addMatrix(a11, a22);
        strassOperations += operationCost;
        double[][] s6 = addMatrix(b11, b22);
        strassOperations += operationCost;
        double[][] s7 = subMatrix(a12, a22);
        strassOperations += operationCost;
        double[][] s8 = addMatrix(b21, b22);
        strassOperations += operationCost;
        double[][] s9 = subMatrix(a11, a21);
        strassOperations += operationCost;
        double[][] s10 = addMatrix(b11, b12);
        strassOperations += operationCost;

       // Μέρος 3
        double[][] p1 = strassen(a11, s1);
        double[][] p2 = strassen(s2, b22);
        double[][] p3 = strassen(s3, b11);
        double[][] p4 = strassen(a22, s4);
        double[][] p5 = strassen(s5, s6);
        double[][] p6 = strassen(s7, s8);
        double[][] p7 = strassen(s9, s10);

        operationCost = p1.length * p1.length;

        // Μέρος 4
        double[][] c11 = addMatrix(subMatrix(addMatrix(p5, p4), p2), p6);
        strassOperations += operationCost + operationCost + operationCost;
        double[][] c12 = addMatrix(p1, p2);
        strassOperations += operationCost;
        double[][] c21 = addMatrix(p3, p4);
        strassOperations += operationCost;
        double[][] c22 = subMatrix(subMatrix(addMatrix(p5, p1), p3), p7);
        strassOperations += operationCost + operationCost + operationCost;
        
     // Τελικός πίνακα
        for (int i = 0; i < c11.length; i++) {
            System.arraycopy(c11[i], 0, productMat[i], 0, c11.length);
        }
        for (int i = 0; i < c12.length; i++) {
            System.arraycopy(c12[i], 0, productMat[i], n / 2, c12.length);
        }
        for (int i = 0; i < c21.length; i++) {
            System.arraycopy(c21[i], 0, productMat[i + n / 2], 0, c21.length);
        }
        for (int i = 0; i < c22.length; i++) {
            System.arraycopy(c22[i], 0, productMat[i + n / 2], n / 2, c22.length);
        }
        return productMat;
    }

}