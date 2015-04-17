#include <stdio.h>
#include <wcslib/wcs.h>
#include "opm.h"
#include "data/small_ref.h"

int
main(int argc, char* argv[])
{
  wcsprm* wcs = new wcsprm;

  wcsini(0, 2, wcs);
  wcs->crpix[0] = 512;
  wcs->crpix[1] = 512;
  wcs->cdelt[0] = -0.00026278;
  wcs->cdelt[1] =  0.00026278;
  wcs->crval[0] =  160.;
  wcs->crval[1] =   35.;
  sprintf(wcs->cunit[0],"deg");
  sprintf(wcs->cunit[1],"deg");
  sprintf(wcs->ctype[0],"RA---TAN");
  sprintf(wcs->ctype[1],"DEC--TAN");
  sprintf(wcs->wcsname,"FK5");
  wcsset(wcs);
  wcsprt(wcs);

  double lon = 160.0;
  double lat =  35.0;
  for (int j=0; j<5; j++)
    for (int i=0; i<5; i++) {
      int status;
      double phi, theta;
      double pixcrd[2], imgcrd[2];
      double a = lon + 0.1*i - 2.0;
      double d = lat + 0.1*i - 2.0;
      double world[2] = {a, d};
      wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &status);
      printf("(%6.2lf,%6.2lf)\n", a, d);
      printf("    => (%8.2lf,%8.2lf) native\n", phi, theta);
      printf("    => (%8.2lf,%8.2lf) image\n", imgcrd[0], imgcrd[1]);
      printf("    => (%8.1lf,%8.1lf) pixel\n", pixcrd[0], pixcrd[1]);
    }

  printf("\n\n");
  auto &ref = opm_test::small_ref;
  for (int i=0; i<10; i++) {
      int status;
      double phi, theta;
      double pixcrd[2], imgcrd[2];
      double world[2] = {ref[i].x, ref[i].y};
      wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &status);
      printf("(%12.8lf,%12.8lf)\n", world[0], world[1]);
      printf("    => (%16.8lf,%16.8lf) native\n", phi, theta);
      printf("    => (%16.8lf,%16.8lf) native\n",
             cos(phi*M_PI/180.)*cos(theta*M_PI/180.)/M_PI*180.,
             sin(phi*M_PI/180.)*cos(theta*M_PI/180.)/M_PI*180.);
      printf("    => (%16.8lf,%16.8lf) image\n", imgcrd[0], imgcrd[1]);
      printf("    => (%16.8lf,%16.8lf) logical\n", 
             imgcrd[0]/wcs->cdelt[0]+wcs->crpix[0],
             imgcrd[1]/wcs->cdelt[1]+wcs->crpix[1]);
      printf("    => (%16.8lf,%16.8lf) pixel\n", pixcrd[0], pixcrd[1]);
  }

  delete wcs;
  return 0;
}
