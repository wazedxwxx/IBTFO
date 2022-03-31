#include "Level_Set.H"
#define Index(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
#include <iostream>
#include <math.h>
#include "PngProcess.H"
#include "ParamReader.H"
using namespace std;
// 0 x_coord 1 y_coord 2 phi 3 n_x 4 n_y 5 cell type

void Level_Set(char *filename,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD)
{

    int image_width;
    int image_height;
    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    const double lo_x = para.get("lo_x", 0);   // x direction grid number in computational domain
    const double lo_y = para.get("lo_y", 0);   // y direction grid number in computational domain
    const double hi_x = para.get("hi_x", 1);   // x direction grid number in computational domain
    const double hi_y = para.get("hi_y", 0.4); // y direction grid number in computational domain

    const double Psy_L = hi_x - lo_x;
    const double Psy_H = hi_y - lo_y;

    cout << " ====  read image ====" << endl;

    png_bytep *row_pointers = read_png_file("/home/ubuntu/IB2TF/Test/PNGGeo/A.png", &image_width, &image_height);

    double *TempImage = new double[image_width * image_height * 6];
// Caculate pixel value
#pragma acc parallel loop
    for (int y = 0; y < image_height; y++)
    {
        png_bytep row = row_pointers[y];
#pragma acc loop
        for (int x = 0; x < image_width; x++)
        {
            png_bytep px = &(row[x * 4]);
            double R_pixel = (double)px[0];
            double G_pixel = (double)px[1];
            double B_pixel = (double)px[2];
            TempImage[Index(x, image_height - y - 1, 0, image_width)] = ((R_pixel + G_pixel + B_pixel) > 0) * 2 - 1;
        }
    }

    // Caculate pixel gradient
    int edg_num = 0;

#pragma acc parallel loop 
    for (int y = 1; y < image_height - 1; y++)
    {
#pragma acc loop
        for (int x = 1; x < image_width - 1; x++)
        {
            double dx = 0;
            double dy = 0;
            if (x > 0)
                dx = TempImage[Index(x, y, 0, image_width)] - TempImage[Index(x - 1, y, 0, image_width)];
            else
                dx = TempImage[Index(x + 1, y, 0, image_width)] - TempImage[Index(x, y, 0, image_width)];

            if (y > 0)
                dy = TempImage[Index(x, y, 0, image_width)] - TempImage[Index(x, y - 1, 0, image_width)];
            else
                dy = TempImage[Index(x, y + 1, 0, image_width)] - TempImage[Index(x, y, 0, image_width)];
            TempImage[Index(x, y, 1, image_width)] = (abs(dx) > abs(dy)) ? abs(dx) : abs(dy);
            if (TempImage[Index(x, y, 1, image_width)] > 0.1){
            #pragma acc atomic update
                edg_num++;
            }

        }
    }

    double *edg_point = new double[edg_num * 2];
    int edg_point_idx = 0;

    for (int y = 0; y < image_height; y++)
    {
        for (int x = 0; x < image_width; x++)
        {
            double dx = (Psy_L / N_x) * (N_x + 2 * num_ghost_cell) / image_width;
            double dy = (Psy_H / N_y) * (N_y + 2 * num_ghost_cell) / image_height;
            if (TempImage[Index(x, y, 1, image_width)] > 0.1)
            {
                edg_point[edg_point_idx] = -num_ghost_cell * (Psy_L / N_x) + dx * x;
                edg_point[edg_num + edg_point_idx] = -num_ghost_cell * (Psy_H / N_y) + dy * y;
                edg_point_idx++; 
            }
        }
    }


#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            int x = (i + 0.5) * (image_width) / (N_x + 2 * num_ghost_cell);
            int y = (j + 0.5) * (image_height) / (N_y + 2 * num_ghost_cell);
            // cout << i << " " << j << " " << TempImage[Index(x, y, 0, image_width)] << endl;
            XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = TempImage[Index(x, y, 0, image_width)]; // Caculate cell type
        }
    }

    double min_distance;
#pragma acc parallel loop reduction(min \
                                    : min_distance)
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            /*            XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)];
                        int count_factor = 0;
            #pragma acc loop
                        for (int k1 = -num_ghost_cell - 1; k1 < num_ghost_cell + 2; k1++)
                        {
            #pragma acc loop
                            for (int k2 = -num_ghost_cell - 1; k2 < num_ghost_cell + 2; k2++)
                            {
                                count_factor = count_factor + XYCOORD[Index(i + k1, j + k2, 2, N_x + 2 * num_ghost_cell)];
                            }
                        }

                        if (count_factor != (2 * num_ghost_cell + 3) * (2 * num_ghost_cell + 3) * XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)])
                        {*/
            min_distance = 100;
            double distance;
            int pn_factor = 1;
#pragma acc loop
            for (int k3 = 0; k3 < edg_num; k3++)
            {
                distance = sqrt((edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                    (edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                    (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]));

                if (distance < min_distance)
                    min_distance = distance;
            }
            pn_factor = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] >= 0 ? 1 : -1;
            XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = pn_factor * min_distance;
            //}
        }
    }

    /*#pragma acc parallel loop
        for (int i = 0; i < num_ghost_cell + 1; i++)
        {
    #pragma acc loop
            for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
            {
                min_distance = 1000;
                for (int k3 = 0; k3 < edg_num; k3++)
                {
                    double distance = sqrt((edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                           (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]));

                    if (distance < min_distance)
                        min_distance = distance;
                }
                int pn_factor = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] > 0 ? 1 : -1;
                XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = pn_factor * min_distance;
            }
        }

    #pragma acc parallel loop
        for (int i = N_x + num_ghost_cell - 1; i < N_x + 2 * num_ghost_cell; i++)
        {
    #pragma acc loop
            for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
            {
                min_distance = 1000;
                for (int k3 = 0; k3 < edg_num; k3++)
                {
                    double distance = sqrt((edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                           (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]));

                    if (distance < min_distance)
                        min_distance = distance;
                }
                int pn_factor = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] > 0 ? 1 : -1;
                XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = pn_factor * min_distance;
            }
        }

    #pragma acc parallel loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
    #pragma acc loop
            for (int j = 0; j < num_ghost_cell + 1; j++)
            {
                min_distance = 1000;
                for (int k3 = 0; k3 < edg_num; k3++)
                {
                    double distance = sqrt((edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                           (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]));

                    if (distance < min_distance)
                        min_distance = distance;
                }
                int pn_factor = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] > 0 ? 1 : -1;
                XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = pn_factor * min_distance;
            }
        }

    #pragma acc parallel loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
    #pragma acc loop
            for (int j = N_y + num_ghost_cell - 1; j < N_x + 2 * num_ghost_cell; j++)
            {
                min_distance = 1000;
                for (int k3 = 0; k3 < edg_num; k3++)
                {
                    double distance = sqrt((edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3] - XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                           (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                               (edg_point[k3 + edg_num] - XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]));

                    if (distance < min_distance)
                        min_distance = distance;
                }
                int pn_factor = XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] > 0 ? 1 : -1;
                XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = pn_factor * min_distance;
            }
        }

    #pragma acc parallel loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
    #pragma acc loop
            for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
            {
                XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)];
                //cout << i << " " << j << " " << XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] << " " << XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] << " " << XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] << endl;
            }
        }*/

#pragma acc parallel loop
    for (int i = 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 1; j < N_y + 2 * num_ghost_cell - 1; j++)
        {
            XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = (XYCOORD[Index(i + 1, j, 2, N_x + 2 * num_ghost_cell)] -
                                                                 XYCOORD[Index(i - 1, j, 2, N_x + 2 * num_ghost_cell)]) /
                                                                ((XYCOORD[Index(i + 1, j, 0, N_x + 2 * num_ghost_cell)] -
                                                                  XYCOORD[Index(i - 1, j, 0, N_x + 2 * num_ghost_cell)])); // N_x
            XYCOORD[Index(i, j, 4, N_x + 2 * num_ghost_cell)] = (XYCOORD[Index(i, j + 1, 2, N_x + 2 * num_ghost_cell)] -
                                                                 XYCOORD[Index(i, j - 1, 2, N_x + 2 * num_ghost_cell)]) /
                                                                ((XYCOORD[Index(i, j + 1, 1, N_x + 2 * num_ghost_cell)] -
                                                                  XYCOORD[Index(i, j - 1, 1, N_x + 2 * num_ghost_cell)])); // N_y
        }
    }

    // Identify each cell type, fluid cell(0), solid cell(1), and ghost cell(2)

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] >= 0)
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 0; // fluid cell
            else
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 1; // solid cell
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] == 0)
            {
                for (int k = 1; k < num_ghost_cell + 1; k++)
                {
                    XYCOORD[Index(i + k, j, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i + k, j, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i - k, j, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i - k, j, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i, j + k, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i, j + k, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i, j - k, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i, j - k, 5, N_x + 2 * num_ghost_cell)];
                }
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = N_x + num_ghost_cell; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = N_y + num_ghost_cell; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = N_y + num_ghost_cell; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

    delete[] edg_point;
    delete[] TempImage;
}