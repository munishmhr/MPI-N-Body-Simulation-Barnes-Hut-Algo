#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

typedef struct QuadTree {
    double NWX[2];
    double NWY[2];
    double NEX[2];
    double NEY[2];
    double SWX[2];
    double SWY[2];
    double SEX[2];
    double SEY[2];
    int depth;
    double xAxis;
    double yAxis;
    double centerOfMass;
    double singleBodyMass;
    int isExternal;
    struct QuadTree* northWest;
    struct QuadTree* northEast;
    struct QuadTree* southWest;
    struct QuadTree* southEast;
} node;
typedef node *tree;

int myid, size;
double starttime, endtime;
FILE *fe;

void sendNodesToCreateQuadTree(double *mass, double *x, double *y, double *vX, double *vY, int noOfPoints, double mX, double mY);
void insertNode(tree temp, tree root, double mass, double x, double y);
void rangeReturn(double x1, double x2, double y1, double y2, tree temp);
void depthSet(tree root, tree temp);
int checkRangeNE(tree p, double x, double y, int o);
int checkRangeNW(tree p, double x, double y, int o);
int checkRangeSE(tree p, double x, double y, int o);
int checkRangeSW(tree p, double x, double y, int o);
void treeTraversal(tree root);
void downGradeAndNextCall(tree movingNode, tree downgrade, tree temp, double mass, double x, double y);
double *UpdateExternalNodes(tree root, double totalMass, double xAxis, double yAxis);
void UpdatedTreeTraversal(tree root);
double ForceCalculation(double mass, double x, double y, double vX, double vY, double mX, double mY, tree root);

void rangeReturn(double x1, double x2, double y1, double y2, tree temp) {
    temp->NEX[0] = ((x1)+(x2)) / 2;
    temp->NEX[1] = x2;
    temp->NEY[0] = ((y1)+(y2)) / 2;
    temp->NEY[1] = y2;
    temp->NWX[0] = x1;
    temp->NWX[1] = ((x1)+(x2)) / 2;
    temp->NWY[0] = ((y1)+(y2)) / 2;
    temp->NWY[1] = y2;
    temp->SWX[0] = x1;
    temp->SWX[1] = ((x1)+(x2)) / 2;
    temp->SWY[0] = y1;
    temp->SWY[1] = ((y1)+(y2)) / 2;
    temp->SEX[0] = ((x1)+(x2)) / 2;
    temp->SEX[1] = x2;
    temp->SEY[0] = y1;
    temp->SEY[1] = ((y1)+(y2)) / 2;
}

void depthSet(tree root, tree temp) {
    temp->depth = root->depth + 1;
}

void insertNode(tree temp, tree movingNode, double mass, double x, double y) {
   // printf("\npoints %lf %lf ",x,y);
    if (movingNode->singleBodyMass == 0) {
        /*
         * Code for NorthEast block
         */
        if (checkRangeNE(movingNode, x, y, 1) == 1) {
            if (movingNode->northEast == NULL) {
                movingNode->northEast = temp;
                rangeReturn(movingNode->NEX[0], movingNode->NEX[1], movingNode->NEY[0], movingNode->NEY[1], temp);
                depthSet(movingNode, movingNode->northEast);
                //printf("\nDepth is %d Inside NE %lf %lf %lf %lf %lf %lf", movingNode->northEast->depth, movingNode->NEX[0], movingNode->NEX[1], movingNode->NEY[0], movingNode->NEY[1], x, y);
            } else {
                if (movingNode->northEast->singleBodyMass == 0) {
                    movingNode = movingNode->northEast;
                    insertNode(temp, movingNode, mass, x, y);
                } else if (movingNode->northEast->singleBodyMass != 0) {
                    tree downgrade = (tree) malloc(sizeof (node));
                    downgrade->singleBodyMass = movingNode->northEast->singleBodyMass;
                    downgrade->xAxis = movingNode->northEast->xAxis;
                    downgrade->yAxis = movingNode->northEast->yAxis;
                    movingNode->northEast->singleBodyMass = 0;
                    movingNode->northEast->isExternal = 1;
                    downGradeAndNextCall(movingNode->northEast, downgrade, temp, mass, x, y);
                }
            }
        } else

            /*
             * Code for NorthWest block
             */
            if (checkRangeNW(movingNode, x, y, 1) == 1) {
            if (movingNode->northWest == NULL) {
                movingNode->northWest = temp;
                rangeReturn(movingNode->NWX[0], movingNode->NWX[1], movingNode->NWY[0], movingNode->NWY[1], temp);
                depthSet(movingNode, movingNode->northWest);
                //printf("\nDepth is %d Inside NW %lf %lf %lf %lf %lf %lf", movingNode->northWest->depth, movingNode->NWX[0], movingNode->NWX[1], movingNode->NWY[0], movingNode->NWY[1], x, y);

            } else {
                if (movingNode->northWest->singleBodyMass == 0) {
                    movingNode = movingNode->northWest;
                    insertNode(temp, movingNode, mass, x, y);
                } else if (movingNode->northWest->singleBodyMass != 0) {
                    tree downgrade = (tree) malloc(sizeof (node));
                    downgrade->singleBodyMass = movingNode->northWest->singleBodyMass;
                    downgrade->xAxis = movingNode->northWest->xAxis;
                    downgrade->yAxis = movingNode->northWest->yAxis;
                    movingNode->northWest->singleBodyMass = 0;
                    movingNode->northWest->isExternal = 1;
                    downGradeAndNextCall(movingNode->northWest, downgrade, temp, mass, x, y);

                }

            }

        } else

            /*
             * Code for SouthWest block
             */
            if (checkRangeSW(movingNode, x, y, 1) == 1) {
            if (movingNode->southWest == NULL) {
                movingNode->southWest = temp;
                rangeReturn(movingNode->SWX[0], movingNode->SWX[1], movingNode->SWY[0], movingNode->SWY[1], temp);
                depthSet(movingNode, movingNode->southWest);
                //printf("\nDepth is %d Inside SW %lf %lf %lf %lf %lf %lf", movingNode->southWest->depth, movingNode->SWX[0], movingNode->SWX[1], movingNode->SWY[0], movingNode->SWY[1], x, y);

            } else {
                if (movingNode->southWest->singleBodyMass == 0) {
                    movingNode = movingNode->southWest;
                    insertNode(temp, movingNode, mass, x, y);
                } else if (movingNode->southWest->singleBodyMass != 0) {
                    tree downgrade = (tree) malloc(sizeof (node));
                    downgrade->singleBodyMass = movingNode->southWest->singleBodyMass;
                    downgrade->xAxis = movingNode->southWest->xAxis;
                    downgrade->yAxis = movingNode->southWest->yAxis;
                    movingNode->southWest->singleBodyMass = 0;
                    movingNode->southWest->isExternal = 1;
                    downGradeAndNextCall(movingNode->southWest, downgrade, temp, mass, x, y);
                }
            }
        } else
            /*
             * Code for SsouthEast block
             */
            if (checkRangeSE(movingNode, x, y, 1) == 1) {
            if (movingNode->southEast == NULL) {
                movingNode->southEast = temp;
                rangeReturn(movingNode->SEX[0], movingNode->SEX[1], movingNode->SEY[0], movingNode->SEY[1], temp);
                depthSet(movingNode, movingNode->southEast);
                //printf("\nDepth is %d Inside SE %lf %lf %lf %lf %lf %lf", movingNode->southEast->depth, movingNode->SEX[0], movingNode->SEX[1], movingNode->SEY[0], movingNode->SEY[1], x, y);
            } else {
                if (movingNode->southEast->singleBodyMass == 0) {
                    movingNode = movingNode->southEast;
                    insertNode(temp, movingNode, mass, x, y);
                } else if (movingNode->southEast->singleBodyMass != 0) {
                    tree downgrade = (tree) malloc(sizeof (node));
                    downgrade->singleBodyMass = movingNode->southEast->singleBodyMass;
                    downgrade->xAxis = movingNode->southEast->xAxis;
                    downgrade->yAxis = movingNode->southEast->yAxis;
                    movingNode->southEast->singleBodyMass = 0;
                    movingNode->southEast->isExternal = 1;
                    downGradeAndNextCall(movingNode->southEast, downgrade, temp, mass, x, y);
                }
            }
        }
    }
}

double ForceCalculation(double mass, double x, double y, double vX, double vY, double mX, double mY, tree root) {
    double force = 0;
    if (root != NULL) {
        double r = (((x - root->xAxis)*(x - root->xAxis))+((y - root->yAxis)*(y - root->yAxis)));
        if (root->isExternal == 0 && root->xAxis != x && root->yAxis != y) {
            force = (root->singleBodyMass * mass) / r;
        } else {
            double c = ((2 * mX)*(2 * mX)) + (2 * mY)*(2 * mY) / (root->depth + 1);
            if ((c / r) < .2 && root != NULL) {
                force = (root->singleBodyMass * mass) / r;
                force = .993;
            } else {
                force = force + ForceCalculation(mass, x, y, vX, vY, mX, mY, root->northEast);
                force = force + ForceCalculation(mass, x, y, vX, vY, mX, mY, root->northWest);
                force = force + ForceCalculation(mass, x, y, vX, vY, mX, mY, root->southWest);
                force = force + ForceCalculation(mass, x, y, vX, vY, mX, mY, root->southEast);
            }
        }
    }
    return force;
}

double *UpdateExternalNodes(tree root, double totalMass, double xAxis, double yAxis) {
if(root != NULL){
//printf("\ninside 1 ");
    double *ar = (double*)realloc(ar,3 * sizeof (*ar));

    if (root->northEast != NULL && root->northEast->singleBodyMass != 0) {
       // printf("\ninside 1 %lf",root->northEast->singleBodyMass);
        totalMass = root->northEast->singleBodyMass + totalMass;
    } else {
        if (root->northEast != NULL && root->northEast->singleBodyMass == 0) {
                //printf("\ninside 11 %lf",root->northEast->singleBodyMass);
            //ar = UpdateExternalNodes(root->northEast, 0, 0, 0);
            totalMass = ar[0] + totalMass;
//	if(root->northEast!=NULL){
//		root = root->northEast;

//	}
	            }
    }

    //printf("\nInside 2 ");

    if (root->northWest != NULL && root->northWest->singleBodyMass != 0) {
        //printf("\ninside 2 %lf",root->northWest->singleBodyMass);
        totalMass = root->northWest->singleBodyMass + totalMass;

    } else {
        if (root->northWest != NULL && root->northWest->singleBodyMass == 0) {
               // printf("\ninside 22 %lf",root->northWest->singleBodyMass);
            //ar = UpdateExternalNodes(root->northWest, 0, 0, 0);
            totalMass = ar[0] + totalMass;
//if(root->northWest!=NULL){
//		root = root->northWest;
//}
        }
    }
   //printf("\ninside 3 ");
    if (root->southEast != NULL && root->southEast->singleBodyMass != 0) {
            //printf("\ninside 3 %lf",root->southEast->singleBodyMass);
        totalMass = root->southEast->singleBodyMass + totalMass;
    } else {
        if (root->southEast != NULL && root->southEast->singleBodyMass == 0) {
               // printf("\ninside 33 %lf",root->southEast->singleBodyMass);
           // ar = UpdateExternalNodes(root->southEast, 0, 0, 0);
            totalMass = ar[0] + totalMass;
//if(root->southEast!=NULL){
//root = root->southEast;
//
//}
        }
    }
    //printf("\ninside 4 ");
    if (root->southWest != NULL && root->southWest->singleBodyMass != 0) {
        //printf("\ninside 4 %lf",root->southWest->singleBodyMass);

      totalMass = root->southWest->singleBodyMass + totalMass;

    } else {
        if (root->southWest != NULL && root->southWest->singleBodyMass == 0) {
                //printf("\ninside 44 %lf",root->southWest->singleBodyMass);
             //ar = UpdateExternalNodes(root->southWest, 0, 0, 0);
             totalMass = ar[0] +totalMass;
//if(root->southWest!=NULL){
//	   root = root->southWest;
//}
        }
    }

    root->singleBodyMass = totalMass;

    ar[0] = totalMass;
    //printf("\n%lf",totalMass);
    return ar;
}
}

int main(int argc, char** argv) {
    /*
        Reading from file and storing mass x and y axis in arrays
     */
    int myid, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0){

        starttime = MPI_Wtime();
    }
        FILE *fp;
        int i = 0;
        double maxY = 0;
        double maxX = 0;
        fp = fopen("tab4000", "r");
        double mass[200000];
        double x[200000];
        double y[200000];
        double z[200000];
        double vX[200000];
        double vY[200000];
        double vZ[200000];

        if (fp != NULL) {
            while (!feof(fp)) {
                fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &mass[i], &x[i], &y[i], &z[i], &vX[i], &vY[i], &vZ[i]);
                if (maxX < labs(x[i])) {
                    maxX = labs(x[i]);
                }
                if (maxY < labs(y[i])) {
                    maxY = labs(y[i]);
                }
                i++;
            }
        }

    int noOfPoints = i - 1;
    double mX = maxX;
    double mY = maxY;
    int f;

    for (f = 1; f <= 1; f++) {
        int count1;
        double force = 0;
        mX = 2 * mX;
        mY = 2 * mY;
        tree root = NULL;
        root = (tree) malloc(sizeof (node));
        tree movingNode = NULL;
        movingNode = (tree) malloc(sizeof (node));
        root->NEX[0] = 0;
        root->NEX[1] = 1 + mX / 2;
        root->NEY[0] = 0;
        root->NEY[1] = 1 + mY / 2;
        root->NWX[0] = -1 - mX / 2;
        root->NWX[1] = 0;
        root->NWY[0] = 0;
        root->NWY[1] = 1 + mY / 2;
        root->SWX[0] = -1 - mX / 2;
        root->SWX[1] = 0;
        root->SWY[0] = -1 - mY / 2;
        root->SWY[1] = 0;
        root->SEX[0] = 0;
        root->SEX[1] = 1 + mX / 2;
        root->SEY[0] = -1 - mY / 2;
        root->SEY[1] = 0;
        root->northEast = NULL;
        root->northWest = NULL;
        root->southEast = NULL;
        root->southWest = NULL;
        root->singleBodyMass = 0;
        root->centerOfMass = 0;
        root->isExternal = 1;

        for (count1 = 0; count1 < noOfPoints; count1++) {
            tree temp = NULL;
            temp = (tree) malloc(sizeof (node));
            temp->northEast = NULL;
            temp->northWest = NULL;
            temp->southEast = NULL;
            temp->southWest = NULL;
            temp->xAxis = x[count1];
            temp->yAxis = y[count1];
            temp->singleBodyMass = mass[count1];
            temp->centerOfMass = 0;
            temp->isExternal = 0;
            movingNode = root;
            insertNode(temp, movingNode, mass[count1], x[count1], y[count1]);
        }
        //printf("\n\n**************** Printing Tree Quadron Wise ****************\n\n");
        //treeTraversal(root);

        //printf("\n\n**************** Updating External Nodes ****************\n\n");
        UpdateExternalNodes(root, 0, 0, 0);
//        printf("\n\n**************** Calculating Force ****************\n\n");
//
        double mXN = 0, mYN = 0;
        double forceArryX[noOfPoints];
        double forceArryY[noOfPoints];
        int subPointsProcessor = noOfPoints/size;
        double subForceArrayX[subPointsProcessor];
        double subForceArrayY[subPointsProcessor];
        int g = 0 ;
        //printf("\ninside loop 1");
        for (count1 = (subPointsProcessor*myid); count1 < (subPointsProcessor*myid)+subPointsProcessor; count1++) {
                //printf("\ninside loop %d %d %d",(subPointsProcessor*myid),(subPointsProcessor*myid)+subPointsProcessor,myid);
            subForceArrayX[g] = ForceCalculation(mass[count1], x[count1], y[count1], vX[count1], vY[count1], mX, mY, root);
            g++;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&subForceArrayX, subPointsProcessor, MPI_DOUBLE, forceArryX, subPointsProcessor, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(forceArryX,noOfPoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
        //MPI_Bcast(forceArryY,noOfPoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (count1 = 0; count1 < noOfPoints; count1++) {
            vX[count1] = forceArryX[count1] / mass[count1] + vX[count1];
            vY[count1] = forceArryX[count1] / mass[count1] + vY[count1];
            x[count1] = vX[count1] + x[count1];
            y[count1] = vY[count1] + y[count1];

            if (mXN < labs(x[count1])) {
                mXN = labs(x[count1]);
            }
            if (mYN < labs(y[count1])) {
                mYN = labs(y[count1]);
            }
        }
        mX = mXN;
        mY = mYN;
    }
       MPI_Barrier(MPI_COMM_WORLD);

    // printf("\ntime : %lf",(double)(endtime));
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0){
        endtime = MPI_Wtime();
//        fe = fopen("output.txt","a+");
//        fprintf(fe,"%d ", size);
        printf("\ntime : %lf",(double)(endtime - starttime));
//        fprintf(fe,"%f", (double)(endtime - starttime));
//        fclose(fe);
    }
    fclose(fp);
     MPI_Finalize();
}

void downGradeAndNextCall(tree movingNode, tree downgrade, tree temp, double mass, double x, double y) {
    if (checkRangeNE(movingNode, downgrade->xAxis, downgrade->yAxis, 0)) {
        movingNode ->northEast = downgrade;
        rangeReturn(movingNode->NEX[0], movingNode->NEX[1], movingNode->NEY[0], movingNode->NEY[1], downgrade);
        depthSet(movingNode, movingNode ->northEast);
        insertNode(temp, movingNode, mass, x, y);
    }
    if (checkRangeNW(movingNode, downgrade->xAxis, downgrade->yAxis, 0)) {
        movingNode ->northWest = downgrade;
        rangeReturn(movingNode->NWX[0], movingNode->NWX[1], movingNode->NWY[0], movingNode->NWY[1], downgrade);
        depthSet(movingNode, movingNode->northWest);
        insertNode(temp, movingNode, mass, x, y);
    }
    if (checkRangeSW(movingNode, downgrade->xAxis, downgrade->yAxis, 0)) {
        movingNode ->southWest = downgrade;
        rangeReturn(movingNode->SWX[0], movingNode->SWX[1], movingNode->SWY[0], movingNode->SWY[1], downgrade);
        depthSet(movingNode, movingNode ->southWest);
        insertNode(temp, movingNode, mass, x, y);
    }
    if (checkRangeSE(movingNode, downgrade->xAxis, downgrade->yAxis, 0)) {
        movingNode ->southEast = downgrade;
        rangeReturn(movingNode->SEX[0], movingNode->SEX[1], movingNode->SEY[0], movingNode->SEY[1], downgrade);
        depthSet(movingNode, movingNode ->southEast);
        insertNode(temp, movingNode, mass, x, y);
    }
}

int checkRangeNE(tree mN, double x, double y, int o) {
    if (mN->NEX[0] <= x && x < mN->NEX[1] && mN->NEY[0] <= y && y < mN->NEY[1]) {
        return 1;
    }
    return 0;
}

int checkRangeNW(tree mN, double x, double y, int o) {
    if (mN->NWX[0] <= x && x < mN->NWX[1] && mN->NWY[0] <= y && y < mN->NWY[1]) {
        return 1;
    }
    return 0;
}

int checkRangeSW(tree mN, double x, double y, int o) {
    if (mN->SWX[0] <= x && x < mN->SWX[1] && mN->SWY[0] <= y && y < mN->SWY[1]) {
        return 1;
    }
    return 0;
}

int checkRangeSE(tree mN, double x, double y, int o) {
    if (mN->SEX[0] <= x && x < mN->SEX[1] && mN->SEY[0] <= y && y < mN->SEY[1]) {
        return 1;
    }
    return 0;
}

void treeTraversal(tree root) {
    if (root->northEast != NULL && root->northEast->singleBodyMass != 0) {
        //printf("\nX:%lf,Y:%lf",root->northEast->xAxis,root->northEast->yAxis);
        printf("\nInserted at NE depth: %d,mass:%lf,X range is %lf <= %lf < %lf,Y range is %lf <= %lf < %lf", root->northEast->depth,
                root->northEast->singleBodyMass, root->NEX[0], root->northEast->xAxis, root->NEX[1], root->NEY[0], root->northEast->yAxis, root->NEY[1]);
    } else {
        if (root->northEast != NULL && root->northEast->singleBodyMass == 0) {
            treeTraversal(root->northEast);
        }
    }
    if (root->northWest != NULL && root->northWest->singleBodyMass != 0) {
        // printf("\nX:%lf,Y:%lf",root->northWest->xAxis,root->northWest->yAxis);
        printf("\nInserted at NW depth: %d,mass:%lf,X range is %lf <= %lf < %lf,Y range is %lf <= %lf < %lf", root->northWest->depth,
                root->northWest->singleBodyMass, root->NWX[0], root->northWest->xAxis, root->NWX[1], root->NWY[0], root->northWest->yAxis, root->NWY[1]);
    } else {
        if (root->northWest != NULL && root->northWest->singleBodyMass == 0) {
            treeTraversal(root->northWest);
        }
    }
    if (root->southWest != NULL && root->southWest->singleBodyMass != 0) {
        // printf("\nX:%lf,Y:%lf",root->southWest->xAxis,root->southWest->yAxis);
        printf("\nInserted at SW depth: %d,mass:%lf,X range is %lf <= %lf < %lf,Y range is %lf <= %lf < %lf", root->southWest->depth,
                root->southWest->singleBodyMass, root->SWX[0], root->southWest->xAxis, root->SWX[1], root->SWY[0], root->southWest->yAxis, root->SWY[1]);
    } else {
        if (root->southWest != NULL && root->southWest->singleBodyMass == 0) {
            treeTraversal(root->southWest);
        }
    }
    if (root->southEast != NULL && root->southEast->singleBodyMass != 0) {
        //printf("\nX:%lf,Y:%lf",root->southEast->xAxis,root->southEast->yAxis);
        printf("\nInserted at SE depth: %d,mass:%lf,X range is %lf <= %lf < %lf,Y range is %lf <= %lf < %lf", root->southEast->depth,
                root->southEast->singleBodyMass, root->SEX[0], root->southEast->xAxis, root->SEX[1], root->SEY[0], root->southEast->yAxis, root->SEY[1]);
    } else {
        if (root->southEast != NULL && root->southEast->singleBodyMass == 0) {
            treeTraversal(root->southEast);
        }
    }

}


