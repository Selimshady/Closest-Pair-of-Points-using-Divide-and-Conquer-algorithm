#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SIZE 255


void swap(int*, int*,int); // Change places of x and y values // flag is controlling x and y values
void quickSort(int**,int,int,int); // To sort O(nLogn)
int partition(int**,int,int,int); // Divide the array
int readFileIntoMat(int***); // It returns the size of matrix//And allocate data matrix //gets values //parameter is sending the address of main matrix
void printMatrix(int**,int); // Parameters are matrix and length, prints points to screen
float closestPair(int**,int,int); // Find the closest points between left and right point and returns distance between them, parameters matrix,left and right point.(part)
float compute(int**,int,int); // Find all lengths betweent left and right points // parameters are mat,left and right point(part).
void control(int**,int,int,int,float*); // Last control for two parts of median // parameters matrix,left,median,right,and address of absMin
                                        //So if it needs to be changed, control function changes it.

int main(int argc, char **argv)
{
    int i;
    float min; // To keep min distance
    int **mat; // To keep points
    printf("sample.txt is being read...\n");
	int counts = readFileIntoMat(&mat);//To read file and get the coordinates of points. // and get the size of matrix // send address of matrix
    printf("sample.txt was read.\n");


    
    printf("\n---Unsorted Points---\n\n");
    printMatrix(mat,counts);

    quickSort(mat,0,counts-1,0);// 0 means it should be sorted by x values
    
    printf("\n------Sorted------\n\n");
    printMatrix(mat,counts);
    
    
    printf("\nMin distance is being calculated...\n-------------------------------\n");
    min = closestPair(mat,0,counts-1);//Start the process
    printf("\nMinimum distance is: %lf\n",min);
    
    
    
    
    for(i=0;i<counts;i++)
        free(mat[i]);
    free(mat);
	return 0;
}


float closestPair(int** mat,int l,int r) {
    float leftMin,rightMin,absMin; // to keep min distances // left and right distances
    int median; // Median value
    int n = r-l+1; // size of part of matrix
    if(n> 3) { // check if the size bigger than three
        median = (r+l)/2; //Find median
        leftMin = closestPair(mat,l,median); // calculate the min dist. of left part of median included median
        printf("\n");
        rightMin = closestPair(mat,median+1,r); // calculate the min dist. of right parf of the median
        printf("\n");
        absMin = (leftMin<rightMin) ? leftMin : rightMin; // check which is lesser
        control(mat,l,median,r,&absMin); // control if there is lesser length in 2d distance // if there is, assign it to min 
        printf("\nLeft index: %d,Medium Index: %d,Right index: %d\nLeft min: %f\nRight min: %f\nAbs Min: %f\n--------------------------\n",l,median,r,leftMin,rightMin,absMin);
        
        return absMin; // return the min
    }
    else 
        return compute(mat,l,r); // To find all distances between points
}

void control(int** mat,int l,int med,int r,float* min) {
    int i = med-1; // It keeps the index of left of median
    int j = med+1; // It keeps the index of right of median
    int lasti,lastj;
    int **closePart; // To keep values in distance
    int k,stop; //Stop is important. Stop keeps how many points in the left and right. So It helps us to not calculate again and again for same values
    float dist; // To store distance between two points
    int size = 0; // sizeOf  closePart matrix
    
    lasti= -1;
    lastj= -1;
    
    closePart = (int**)malloc((r-l+1)*sizeof(int*));
    for(k=0;k<=r-l;k++)
        closePart[k] = (int*)malloc(2*sizeof(int));
    
    
   
    while(i>=l && mat[med][0] - mat[i][0] < *min) { // Find the furthest point on left in min distance
        closePart[size][0] = mat[i][0]; // If it is in min distance store it to matrix
        closePart[size][1] = mat[i][1];
        size++; //Index for closePart
        i--; //Index of left part
    }
    
    closePart[size][0] = mat[med][0]; // Store median to where left values finished
    closePart[size][1] = mat[med][1];
    size++;
    stop = size; //The last point of left
    
    
    while(j<=r && mat[j][0] - mat[med][0] < *min) { // Find the furthest point on right in min distance
        closePart[size][0] = mat[j][0]; // If it is in min distance store it to matrix
        closePart[size][1] = mat[j][1];
        size++; //size of close part
        j++; //Index of right part
    }
        
    
    
    if(size>1 && size>stop) { //Checks if there are points in the distance in the right part of the median
        quickSort(closePart,0,size-1,1);// rightest value is flag that shows what based it should be sorted
        for(i=0;i<size-1;i++) { //Index for left part
            j=i+1;
            while(j<size && closePart[j][1] - closePart[i][1] < *min) {//Check next value if it is under min dist.
                    dist = pow(closePart[i][0] - closePart[j][0],2) + pow(closePart[i][1] - closePart[j][1],2); //Compute the distance
                    if(dist < pow(*min,2)) {//If it is lesser than other point, assign it to min
                        lasti = i;
                        lastj = j;
                        *min = sqrt(dist);
                    }
                    j++;
            }
        }
    }
    if(lasti >=0) 
         printf("P[%d][%d] , P[%d][%d]  = %f\n\n",closePart[lasti][0],closePart[lasti][1],closePart[lastj][0],closePart[lastj][1],(*min));
    for(k=0;k<=r-l;k++)
        free(closePart[k]);
    free(closePart);
}

float compute(int** mat,int l,int r) {
    int i,j; // Helps calculate distances between 2 to 3 points
    float min = INFINITY; // Init first distance
    int dist; // To store distance between two points
    for(i=l;i<r;i++) {
        for(j=i+1;j<=r;j++) {
            dist = pow(mat[i][0] - mat[j][0],2) + pow(mat[i][1] - mat[j][1],2); //Compute the distance
            printf("P[%d][%d] , P[%d][%d]  = %f\n\n",mat[i][0],mat[i][1],mat[j][0],mat[j][1],sqrt(dist));
            if(dist < min)//If it is lesser than other point, assign it to min
                min = dist;
        }
    }
    return sqrt(min);//Take square root of it and return the value
}


int readFileIntoMat(int*** mat) {
    FILE *fp; // To keep file pointer
    int tmp[SIZE][2];//Temporary matrix to get length and data
    int size = 0; // Size of points
    int i; // index of matrix
    char* line = (char*)malloc(sizeof(char)*255);//To get the line and values//Allocating memory
    int **tmpMatrix; // Temporary matrix
 
    if((fp = fopen("sample.txt","r")) == NULL) {//Opening the file to read
        printf("FILE ACCSESS ERROR!");
        exit(-1);//If there is a problem, exit the program
    } 
    
    while(fgets(line,SIZE,fp) != NULL) {//Read the line
        sscanf(line, "%d %d" ,&tmp[size][0],&tmp[size][1]); // Parsing two values to matrix
        size++; //Counts the lines
    }
    
    tmpMatrix = (int**)malloc(sizeof(int*)*size); // Allocating the real matrix
    for(i=0;i<size;i++) 
        tmpMatrix[i] = (int*)malloc(sizeof(int)*2);

    for(i=0;i<size;i++) {
        tmpMatrix[i][0] = tmp[i][0];
        tmpMatrix[i][1] = tmp[i][1];
    }
    *mat = tmpMatrix; // To pass address to main matrix
    
    
    fclose(fp); //close the file
    free(line); // free the line address
    return size; //Return the length
}


void printMatrix(int** mat,int counts) {
    int i; //index of matrix 
    for(i=0;i<counts;i++) {
        printf("\t%d %d\n",mat[i][0],mat[i][1]);
    }
}

void quickSort(int** arr,int low,int high,int flag)
{
    int pivot;
    if (low < high)
    {
        /* pivot is partitioning index, arr[pivot] is now
           at right place */
        pivot = partition(arr, low, high,flag);

        quickSort(arr, low, pivot - 1,flag);  // Before pivot //Left part of pivot
        quickSort(arr, pivot + 1, high,flag); // After pivot  //Right part of pivot
    }
}

int partition (int** arr, int low, int high,int flag) 
{ 
    int pivot = arr[high][flag]; // pivot 
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
    int j;
    for (j = low; j <= high - 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j][flag] < pivot) 
        { 
            i++; // increment index of smaller element 
            swap(&arr[i][flag], &arr[j][flag],flag); 
        } 
    } 
    swap(&arr[i + 1][flag], &arr[high][flag],flag); 
    return (i + 1); 
} 

void swap(int* a, int* b,int flag) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
    if(flag == 0) {
        t = *(a+1);
        *(a+1) = *(b+1);
        *(b+1) = t;
    }
    else {
        t = *(a-1);
        *(a-1) = *(b-1);
        *(b-1) = t;
    }
        
}