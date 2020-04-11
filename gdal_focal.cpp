#include <iostream>
#include <cstdint>
#include <clocale>
#include <vector>
#include <gdal_priv.h>
#include <cpl_string.h>

using namespace std;

typedef vector<float> Array;
typedef vector<Array> Matrix;

Matrix getGaussian(int height, int width, double sigma)
{
    Matrix kernel(height, Array(width));
    double sum=0.0;
    int i,j;

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] = exp(-(i*i+j*j)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel[i][j];
        }
    }

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

int main(int argc, char **argv) {

	setlocale(LC_ALL, "Portuguese");

	const char *input_filename = argv[1];
	const char *output_filename = argv[2];

	GDALDataset *in_ds;
	GDALDataset *out_ds;
	GDALDriver *geotiff;
	GDALAllRegister();

	in_ds = (GDALDataset *)GDALOpen(input_filename, GA_ReadOnly);

	int nrows, ncols, nbands;
	float nodata;
	double transform[6];
	const char *proj;
	char **papszOptions = NULL;

	nrows = in_ds->GetRasterBand(1)->GetYSize();
	ncols = in_ds->GetRasterBand(1)->GetXSize();
	nbands = in_ds->GetRasterCount();
	nodata = in_ds->GetRasterBand(1)->GetNoDataValue();
	in_ds->GetGeoTransform(transform);
	proj = in_ds->GetProjectionRef();

	// Matrix filter(3, Array(3));
	Matrix filter = getGaussian(3, 3, 10.0);

	int filterHeight = filter.size();
	int filterWidth = filter[0].size();
	int newImageHeight = nrows - filterHeight + 1;
	int newImageWidth = ncols - filterWidth + 1;
	int d, i, j, h, w;

	geotiff = GetGDALDriverManager()->GetDriverByName("GTiff");
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
	out_ds = geotiff->Create(output_filename, newImageWidth, newImageHeight, 1, GDT_Float32, papszOptions);
	out_ds->SetGeoTransform(transform);
	out_ds->GetRasterBand(1)->SetNoDataValue(nodata);
	out_ds->SetProjection(proj);

	float *input_pixel = (float *)CPLMalloc(sizeof(float));
	float *output_pixel = (float *)CPLMalloc(sizeof(float));

	for (i = 0; i < newImageHeight; i++) {
		for (j = 0; j < newImageWidth; j++) {
			*output_pixel = 0;
			for (h = i; h < i + filterHeight; h++) {
				for (w = j; w < j + filterWidth; w++) {
					in_ds->GetRasterBand(1)->RasterIO(GF_Read, w, h, 1, 1, input_pixel, 1, 1, GDT_Float32, 0, 0);
					*output_pixel += filter[h - i][w - j] + *input_pixel;
				}
			}
			out_ds->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, output_pixel, 1, 1, GDT_Float32, 0, 0);
		}
	}

	CPLFree(input_pixel);
	CPLFree(output_pixel);
	GDALClose(in_ds);
	GDALClose(out_ds);
	GDALDestroyDriverManager();

	cout << "gdal_focal: Process completed on: " << argv[1] << endl;
	cout << "gdal_focal: Your output image " << argv[2] << " is ready!" << endl;
	return 0;
}