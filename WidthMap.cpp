// ==============================================================
// Upon having the skeleton and surface representation of a given
// mitochondria provided by MitoGraph, this routine may be used to
// calculate the the mitochondrial tubules width.
// ==============================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkImageFFT.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>

//#define DEBUG

// Mapping tubules width over the mitochondrial skeleton
void MapOverSkeleton(vtkPolyData *Skeleton, vtkPolyData *Surface, const char FileName[], double *scale);

// Mapping tubules width over the mitochondrial surface
void MapOverSurface(vtkPolyData *Skeleton, vtkPolyData *Surface, const char FileName[], double *scale);

// MitoGraph creates a file called mitograph.files that list
// all the files that were analyzed. This routine uses the same
// file to know what files you are interested in calculating the
// tubules width.
int CalculateTubulesWidth(const char FileName[], double *scale);

// Save PolyData
void SavePolyData(vtkPolyData *PolyData, const char FileName[]);

/* ================================================================
   I/O ROUTINES
=================================================================*/

void SavePolyData(vtkPolyData *PolyData, const char FileName[]) {

    #ifdef DEBUG
        printf("Writing PolyData file %s...\n",FileName);
    #endif

    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer -> SetFileType(VTK_BINARY);
    Writer -> SetInputData(PolyData);
    Writer -> SetFileName(FileName);
    Writer -> Write();
}

/* ================================================================
   -> SURFACE
=================================================================*/

void MapOverSkeleton(vtkPolyData *Skeleton, vtkPolyData *Surface, const char FileName[], double *scale, double *avg_width, double *std_width) {

    #ifdef DEBUG
        printf("Mapping over the skeleton...\n");
    #endif

    vtkIdType N = Skeleton -> GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> Width = vtkSmartPointer<vtkDoubleArray>::New();
    Width -> SetNumberOfComponents(1);
    Width -> SetNumberOfTuples(N);

    #ifdef DEBUG
        printf("Generating point locator...\n");
    #endif

    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree -> SetDataSet(Surface);
    Tree -> BuildLocator();

    int k, n = 3;
    vtkIdType id, idk;
    double r[3], rk[3], d, w, w2;
    vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();
    
    w = w2 = 0.0;
    for (id = 0; id < N; id++) {
        d = 0;
        Skeleton -> GetPoint(id,r);
        Tree -> FindClosestNPoints(n,r,List);
        for (k = 0; k < n; k++) {
            idk = List -> GetId(k);
            Surface -> GetPoint(idk,rk);
            d += sqrt( pow(scale[0]*r[0]-scale[0]*rk[0],2) + pow(scale[0]*r[1]-scale[0]*rk[1],2) + pow(scale[1]*r[2]-scale[1]*rk[2],2) );
        }
        d /= (double)n;
        w += 2*d;
        w2 += pow(2*d,2);
        Width -> SetTuple1(id,2*d);

    }
    Width -> Modified();
    Skeleton -> GetPointData() -> SetScalars(Width);
    *avg_width = w / N;
    *std_width = sqrt(w2-pow(w/N,2)) / N;
}

/* ================================================================
   -> SURFACE
=================================================================*/

void MapOverSurface(vtkPolyData *Skeleton, vtkPolyData *Surface, const char FileName[], double *scale) {

    #ifdef DEBUG
        printf("Mapping over the surface...\n");
    #endif

    vtkIdType N = Surface -> GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> Width = vtkSmartPointer<vtkDoubleArray>::New();
    Width -> SetNumberOfComponents(1);
    Width -> SetNumberOfTuples(N);

    #ifdef DEBUG
        printf("Generating point locator...\n");
    #endif

    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree -> SetDataSet(Skeleton);
    Tree -> BuildLocator();

    int k, n = 3;
    vtkIdType id, idk;
    double r[3], rk[3], d;
    vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();

    for (id = 0; id < N; id++) {
        d = 0;
        Surface -> GetPoint(id,r);
        Tree -> FindClosestNPoints(n,r,List);
        for (k = 0; k < n; k++) {
            idk = List -> GetId(k);
            Skeleton -> GetPoint(idk,rk);
            d += sqrt( pow(scale[0]*r[0]-scale[0]*rk[0],2) + pow(scale[0]*r[1]-scale[0]*rk[1],2) + pow(scale[1]*r[2]-scale[1]*rk[2],2) );
        }
        d /= (double)n;
        Width -> SetTuple1(id,d);
    }
    Width -> Modified();
    Surface -> GetPointData() -> SetScalars(Width);

}

/* ================================================================
   TUBULES WIDTH
=================================================================*/

double CalculateTubulesWidth(const char FileName[], double *scale, double *avg_width, double *std_width) {

    #ifdef DEBUG
        printf("Loading PolyData 1...\n");
    #endif

    // Loading mitochondrial skeleton
    char _fullpath[256];
    sprintf(_fullpath,"%s_skeleton.vtk",FileName);
    vtkSmartPointer<vtkPolyDataReader> SkeletonReader = vtkSmartPointer<vtkPolyDataReader>::New();
    SkeletonReader -> SetFileName(_fullpath);
    int errlog = SkeletonReader -> IsFilePolyData();
    // File cannot be opened
    if (!errlog) {
        printf("File %s does not exist or cannot be opened.\n",_fullpath);
        return -1;
    }
    SkeletonReader -> Update();

    #ifdef DEBUG
        printf("Loading PolyData 2...\n");
    #endif

    // Loading mitochondrial surface
    sprintf(_fullpath,"%s_surface.vtk",FileName);
    vtkSmartPointer<vtkPolyDataReader> SurfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();
    SurfaceReader -> SetFileName(_fullpath);
    errlog = SurfaceReader -> IsFilePolyData();
    // File cannot be opened
    if (!errlog) {
        printf("File %s does not exist or cannot be opened.\n",_fullpath);
        return -1;
    }
    SurfaceReader -> Update();

    vtkPolyData *Skeleton = SkeletonReader -> GetOutput();
    vtkPolyData *Surface  =  SurfaceReader -> GetOutput();

    //MapOverSurface(Skeleton,Surface,FileName,scale);
    MapOverSkeleton(Skeleton,Surface,FileName,scale,avg_width,std_width);

    sprintf(_fullpath,"%s_skeleton-color.vtk",FileName);
    SavePolyData(Skeleton,_fullpath);

    //sprintf(_fullpath,"%s_surface-color.vtk",FileName);
    //SavePolyData(Surface,_fullpath);

    return 0;
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     

    int i;
    char _impath[128];
    sprintf(_impath,"");
    double scale[2] = {1.00,1.00};
    
    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-xy")) {
            scale[0] = (double)atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-z")) {
            scale[1] = (double)atof(argv[i+1]);
        }
    }

    #ifdef DEBUG
        printf("WidthMap...\n");
        printf("File: %s\n",_impath);
        printf("xy = %1.3f, z = %1.3f\n",scale[0],scale[1]);
    #endif

    // Generating list of files to run
    char _cmd[256];
    sprintf(_cmd,"ls %s*_skeleton.vtk | sed -e 's/_skeleton.vtk//' > %swidthmap.files",_impath,_impath);
    system(_cmd);

    char _tifffilename[256];
    char _tifflistpath[128];
    double avg_width, std_width;
    sprintf(_tifflistpath,"%swidthmap.files",_impath);
    FILE *f = fopen(_tifflistpath,"r");
    while (fgets(_tifffilename,256, f) != NULL) {
        _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';
        CalculateTubulesWidth(_tifffilename,scale,&avg_width,&std_width);
        printf("%s\t%1.4f\t%1.4f\n",_tifffilename,avg_width,std_width);
    }
    fclose(f);

}
