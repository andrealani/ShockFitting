// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/LeastSquaresMethod.hh"
#include "MathTools/MovingAverageFilter.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

LeastSquaresMethod::LeastSquaresMethod()
{
}

//--------------------------------------------------------------------------//

LeastSquaresMethod::~LeastSquaresMethod()
{
}

//--------------------------------------------------------------------------//

void LeastSquaresMethod::fitData(vector<double> X,
                                 vector<double> Y,
                                 unsigned polynomialGrad)
{
  cout << "     => LeastSquaresMethod::fitData() -> polynomial order = ";
  cout << polynomialGrad << endl;

  assert(X.size()==Y.size());

  m_nbPoints = X.size();

  bestCoefficients.resize(polynomialGrad+1);

  a.resize(polynomialGrad+1,polynomialGrad+1);
  b.resize(polynomialGrad+1);

  Xnew.resize(m_nbPoints);
  Ynew.resize(m_nbPoints);

  // create object solving the linear system
  Solg<double> SystemSolver;

  // compute the coefficients of the matrix used for the Least Squares Method
  for(unsigned h=0;h<polynomialGrad+1;h++) {
   for(unsigned IPOIN=0;IPOIN<m_nbPoints; IPOIN++) {
    b.at(h) += pow(X.at(IPOIN),h)*Y.at(IPOIN);
   }
   for(unsigned k=0;k<polynomialGrad+1;k++) {
    for(unsigned IPOIN=0;IPOIN<m_nbPoints; IPOIN++) {
     a(h,k) += pow(X.at(IPOIN),h+k);
    }
   }
  }

  // solve the linear system and find the coefficients of the curve which best
  // fit the points distribution
  bestCoefficients = SystemSolver.callSolg(a,b);

  // define the dx used to distribute the X-coordinate
  double dx = (*max_element(X.begin(),X.end())-
               *min_element(X.begin(),X.end()))/
               m_nbPoints;

  // initialize the Y-vector
  Ynew.assign(m_nbPoints,0);

  // compute the new interpolated coordinates
  for(unsigned IPOIN=0;IPOIN<m_nbPoints;IPOIN++) {
   Xnew.at(IPOIN) = *min_element(X.begin(),X.end())+dx*IPOIN;
   for(unsigned j=0;j<polynomialGrad+1;j++) {
    Ynew.at(IPOIN) += bestCoefficients.at(j)*pow(Xnew.at(IPOIN),j); 
   }
  }

  // plot fitted points
  FILE* fittedPoints;
  fittedPoints = fopen("log/LeastSquaresMethod.dat","w");
  
  fprintf(fittedPoints,"%s","TITLE = Fitted Points\n");
  fprintf(fittedPoints,"%s","VARIABLES = \"x0\" \"x1\"\n");
  fprintf(fittedPoints,"%s","ZONE T = \"Fitted points\"\n");
  fprintf(fittedPoints,"%s","STRANDID=0, SOLUTIONTIME=0 ");
  fprintf(fittedPoints,"%s %u %s","I=",m_nbPoints,", J=1, K=1, ZONETYPE=Ordered ");
  fprintf(fittedPoints,"%s","DATAPACKING=POINT\n");
  fprintf(fittedPoints,"%s","DT = (SINGLE, SINGLE)\n");

  for(unsigned IPOIN=0;IPOIN<m_nbPoints;IPOIN++) {
   fprintf(fittedPoints,"%22.14F",Xnew.at(IPOIN));
   fprintf(fittedPoints,"%22.14F",Ynew.at(IPOIN));
   fprintf(fittedPoints,"%s","\n");
  }

  fclose(fittedPoints);
}

//--------------------------------------------------------------------------//

void LeastSquaresMethod::m_fitData(vector<double>& m_X,
                                   vector<double>& m_Y,
                                   unsigned polGrad)
{
  // this method differs from "fitData" because it is used as private
  // method inside the class

  // define object solving linear systems
  Solg<double> SystemSolver;

  vector <double> m_bestCoefficients(polGrad+1);
  Array2D <double> m_a(polGrad+1,polGrad+1);
  vector <double> m_b(polGrad+1,0);

  unsigned nbPoints = m_X.size();

  // compute the coefficients of the matrix used for the Least Squares Method
  for(unsigned h=0;h<polGrad+1;h++) {
   for(unsigned IPOIN=0;IPOIN<nbPoints; IPOIN++) {
    m_b.at(h) += pow(m_X.at(IPOIN),h)*m_Y.at(IPOIN);
   }
   for(unsigned k=0;k<polGrad+1;k++) {
    for(unsigned IPOIN=0;IPOIN<nbPoints; IPOIN++) {
     m_a(h,k) += pow(m_X.at(IPOIN),h+k);
    }
   }
  }

  // solve the linear system and find the coefficients of the curve which best
  // fit the points distribution
  m_bestCoefficients = SystemSolver.callSolg(m_a,m_b);

  // define the dx used to distribute the X-coordinate
  double dx = (*max_element(m_X.begin(),m_X.end())-
               *min_element(m_X.begin(),m_X.end()))/
               nbPoints;

  // make a backup of the X-vector
  vector<double> m_Xold(nbPoints,0);
  for(unsigned IPOIN=0;IPOIN<nbPoints;IPOIN++) {
   m_Xold.at(IPOIN) = m_X.at(IPOIN);
  }

  // initialize the Y-vector
  m_Y.assign(nbPoints,0);

  // compute the new interpolated coordinates
/*  m_X.at(0)=*min_element(m_Xold.begin(),m_Xold.end());

  for(unsigned j=0;j<polGrad+1;j++) {
   m_Y.at(0) +=
    m_bestCoefficients.at(j)*pow(m_X.at(0),j);
  }
  unsigned IPOIN=1;
  while(abs(m_X.at(IPOIN-1) - *max_element(m_Xold.begin(),m_Xold.end()))>=dx) {
   m_X.resize(IPOIN+1);
   m_Y.resize(IPOIN+1);
   // increase the x-coordinates of a dx-step
   m_X.at(IPOIN)=m_X.at(IPOIN-1)+dx;
   for(unsigned j=0;j<polGrad+1;j++) {
    m_Y.at(IPOIN) +=
     m_bestCoefficients.at(j)*pow(m_X.at(IPOIN),j);
   }
   IPOIN++;
  }
*/
  // compute the new interpolated coordinates
  for(unsigned IPOIN=0;IPOIN<nbPoints;IPOIN++) {
   m_X.at(IPOIN)=*min_element(m_Xold.begin(),m_Xold.end())+dx*IPOIN;
   for(unsigned j=0;j<polGrad+1;j++) {
    m_Y.at(IPOIN) +=
     m_bestCoefficients.at(j)*pow(m_X.at(IPOIN),j);
   }
  }

  return;
}

//--------------------------------------------------------------------------//

void LeastSquaresMethod::fitEllipse(vector<double> X,
                                    vector<double> Y)
{
  cout << "     => LeastSquaresMethod::fitEllipse()\n";

  assert(X.size()==Y.size());

  unsigned nbPoints = X.size();
  
  double maxXvalue = *max_element(X.begin(),X.end());
  double minXvalue = *min_element(X.begin(),X.end());

  // number of points used to fit the shock polyline
  unsigned newPoints = 400;

  // compute the deltaX for the x-coordinates
  const double deltaX = (maxXvalue-minXvalue)/(newPoints);

  Xnew.resize(newPoints*2);
  Ynew.resize(newPoints*2);
  bestCoefficients.resize(5);

  a.resize(5,5);
  b.resize(5);

  // create object solving the linear system
  Solg<double> SystemSolver;

  // compute the coefficients of the matrix used for the Least Squares Method
  for(unsigned IPOIN=0;IPOIN<nbPoints;IPOIN++) {
   a(0,0)+= pow(X.at(IPOIN),4);
   a(0,1)+= pow(X.at(IPOIN),3)*Y.at(IPOIN);
   a(0,2)+= pow(X.at(IPOIN),2)*pow(Y.at(IPOIN),2);
   a(0,3)+= pow(X.at(IPOIN),3);
   a(0,4)+= pow(X.at(IPOIN),2)*Y.at(IPOIN);

   a(1,0)+= pow(X.at(IPOIN),3)*Y.at(IPOIN);
   a(1,1)+= pow(X.at(IPOIN),2)*pow(Y.at(IPOIN),2);
   a(1,2)+= X.at(IPOIN)*pow(Y.at(IPOIN),3);
   a(1,3)+= pow(X.at(IPOIN),2)*Y.at(IPOIN);
   a(1,4)+= X.at(IPOIN)*pow(Y.at(IPOIN),2);

   a(2,0)+= pow(X.at(IPOIN),2)*pow(Y.at(IPOIN),2);
   a(2,1)+= X.at(IPOIN)*pow(Y.at(IPOIN),3);
   a(2,2)+= pow(Y.at(IPOIN),4);
   a(2,3)+= X.at(IPOIN)*pow(Y.at(IPOIN),2);
   a(2,4)+= pow(Y.at(IPOIN),3);

   a(3,0)+= pow(X.at(IPOIN),3);
   a(3,1)+= pow(X.at(IPOIN),2)*Y.at(IPOIN);
   a(3,2)+= X.at(IPOIN)*pow(Y.at(IPOIN),2);
   a(3,3)+= pow(X.at(IPOIN),2);
   a(3,4)+= X.at(IPOIN)*Y.at(IPOIN);

   a(4,0)+= pow(X.at(IPOIN),2)*Y.at(IPOIN);
   a(4,1)+= X.at(IPOIN)*pow(Y.at(IPOIN),2);
   a(4,2)+= pow(Y.at(IPOIN),3);
   a(4,3)+= X.at(IPOIN)*Y.at(IPOIN);
   a(4,4)+= pow(Y.at(IPOIN),2);

   b.at(0)+= pow(X.at(IPOIN),2);
   b.at(1)+= X.at(IPOIN)*Y.at(IPOIN);
   b.at(2)+= pow(Y.at(IPOIN),2);
   b.at(3)+= X.at(IPOIN);
   b.at(4)+= Y.at(IPOIN);
  }

  for(unsigned i=0;i<5;i++) { b.at(i)=-b.at(i); }

  // solve the linear system and find the coefficients of the curve which best
  // fit the points distribution
  bestCoefficients = SystemSolver.callSolg(a,b);

  double a = bestCoefficients.at(0);
  double b = bestCoefficients.at(1);
  double c = bestCoefficients.at(2);
  double d = bestCoefficients.at(3);
  double e = bestCoefficients.at(4);

  double x1 = -1.0*(2*b*e-4*c*d) + sqrt((2*b*e-4*c*d)*(2*b*e-4*c*d)-
                                         4.0*(b*b-4*a*c)*(e*e-4*c));
  x1 /= 2*(b*b-4*a*c);
  double x2 = -1.0*(2*b*e-4*c*d) - sqrt((2*b*e-4*c*d)*(2*b*e-4*c*d)-
                                         4.0*(b*b-4*a*c)*(e*e-4*c));
  x2 /= 2*(b*b-4*a*c);

  double m_x;
  m_nbPoints=0;
  for(unsigned IPOIN=0;IPOIN<newPoints;IPOIN++) {
    m_x = minXvalue+deltaX*IPOIN;

    if(m_x>=x1 && m_x<=x2 && m_x<=maxXvalue) {
     Xnew.at(m_nbPoints) = m_x;
     Ynew.at(m_nbPoints) = 
      -1.0*(b*Xnew.at(m_nbPoints)+e)+sqrt((b*Xnew.at(m_nbPoints)+e)*(b*Xnew.at(m_nbPoints)+e)-
                                      4*c*(1+a*Xnew.at(m_nbPoints)*Xnew.at(m_nbPoints)+
                                           d*Xnew.at(m_nbPoints)));
     Ynew.at(m_nbPoints) /= 2.0*c;
     ++m_nbPoints;
    }
    else  {
     cout << "        LeastSquaresMethod:: (!) warning =>";
     cout << " x: " << m_x << " is out of the domain\n";
    }
  }

  unsigned start = m_nbPoints;
  for(unsigned IPOIN=start;IPOIN<newPoints*2;IPOIN++) {
    m_x = Xnew.at(start-1)-deltaX*(m_nbPoints-start);

    if(m_x>=x1 && m_x<=x2 && m_x>=minXvalue) {
     Xnew.at(m_nbPoints) = m_x; 
     Ynew.at(m_nbPoints) = 
      -1.0*(b*Xnew.at(m_nbPoints)+e)-sqrt((b*Xnew.at(m_nbPoints)+e)*(b*Xnew.at(m_nbPoints)+e)-
                                      4*c*(1+a*Xnew.at(m_nbPoints)*Xnew.at(m_nbPoints)+
                                           d*Xnew.at(m_nbPoints)));
     Ynew.at(m_nbPoints) /= 2.0*c;
     ++m_nbPoints;
    }
    else  {
     cout << "        LeastSquaresMethod:: (!) warning =>";
     cout << " x: " << m_x << " is out of the domain\n";
    }
  }

  Xnew.resize(m_nbPoints);
  Ynew.resize(m_nbPoints);
 
  // plot points interpolation 
  FILE* fittedPoints;
  fittedPoints = fopen("log/LeastSquaresMethod.dat","w");

  fprintf(fittedPoints,"%s","TITLE = Fitted Points\n");
  fprintf(fittedPoints,"%s","VARIABLES = \"x0\" \"x1\"\n");
  fprintf(fittedPoints,"%s","ZONE T = \"Fitted points\"\n");
  fprintf(fittedPoints,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
  fprintf(fittedPoints,"%s %u %s","I= ",m_nbPoints,", J=1, K=1, ZONETYPE=Ordered\n");
  fprintf(fittedPoints,"%s","DATAPACKING=POINT\n");
  fprintf(fittedPoints,"%s","DT = (SINGLE, SINGLE)\n");

  for(unsigned IPOIN=0;IPOIN<m_nbPoints;IPOIN++) {
   fprintf(fittedPoints,"%22.14F",Xnew.at(IPOIN));
   fprintf(fittedPoints,"%22.14F",Ynew.at(IPOIN));
   fprintf(fittedPoints,"%s","\n");
  }

  fclose(fittedPoints);
  return;
}

//--------------------------------------------------------------------------//

void LeastSquaresMethod::fitSplittingCurves(vector<double> X,
				       	    vector<double> Y,
                                            vector<unsigned> nbSegments,
                                            vector<unsigned> SegPolynomialOrders,
                                            bool SmoothingOption)
{
  cout << "     => LeastSquaresMethod::fitSplittingCurves()\n";

  assert(X.size()==Y.size());

  m_nbPoints = X.size();

  Xnew.resize(m_nbPoints*2); Ynew.resize(m_nbPoints*2);

  // create object smoothing the curves
  MovingAverageFilter iSmooth;

  // X and Y are sorted from the shock point in bottom of the domain
  // to the shock point on the top of the domain
  const  unsigned nbShockSegmentsX = nbSegments.at(0);
  const  unsigned nbShockSegmentsY = nbSegments.at(1);
  const unsigned nbShockCells = nbShockSegmentsX*nbShockSegmentsY;

  const double minX = *min_element(X.begin(),X.end());
  const double minY = *min_element(Y.begin(),Y.end());

  const double deltaX = abs(*max_element(X.begin(),X.end())-
                            *min_element(X.begin(),X.end()))/
                            nbShockSegmentsX;
  const double deltaY = abs(*max_element(Y.begin(),Y.end())-
                            *min_element(Y.begin(),Y.end()))/
                            nbShockSegmentsY;

  Array2D<unsigned> ptrShCellPoint(m_nbPoints,nbShockCells);
  vector<unsigned> nbCellShockPoints(nbShockCells);
  // define a vector that checks if a point has been already included
  // in a previous cell
  vector<bool> alreadyCounted;
  alreadyCounted.assign(m_nbPoints,false);

  unsigned iShCell=0; unsigned iShCellPoint;

  /// store the points of the distribution according to the cells domains
  /// in which the shock polyline is divided
  for(unsigned ICELL=0;ICELL<nbShockSegmentsX;ICELL++) {
   for(unsigned JCELL=0;JCELL<nbShockSegmentsY;JCELL++) {

    iShCellPoint=0;

    for(unsigned IPOIN=0;IPOIN<m_nbPoints;IPOIN++) {
     if((X.at(IPOIN)-(minX+deltaX*ICELL))     >= -1.e-8 && 
        (X.at(IPOIN)-(minX+deltaX*(ICELL+1))) <= 1.e-8 &&
        (Y.at(IPOIN)-(minY+deltaY*JCELL))     >= -1.e-8 && 
        (Y.at(IPOIN)-(minY+deltaY*(JCELL+1))) <= 1.e-8 ) {
      if(!alreadyCounted.at(IPOIN)) {
       ptrShCellPoint(iShCellPoint,iShCell) = IPOIN;
       ++iShCellPoint;
      }
      alreadyCounted.at(IPOIN) = true;
     }
    }
    if(iShCellPoint!=0) {
     nbCellShockPoints.at(iShCell)=iShCellPoint;
     ++iShCell;
    }
   }
  }

  nbCellShockPoints.resize(iShCell);

  vector<double> m_X;
  vector<double> m_Y; 
  unsigned iShPointNew=0;

  for(unsigned ISHCELL=0;ISHCELL<iShCell;ISHCELL++) {

   m_X.resize(nbCellShockPoints.at(ISHCELL),0);
   m_Y.resize(nbCellShockPoints.at(ISHCELL),0);

   for(unsigned ISHCELLPOINT=0;
                ISHCELLPOINT<nbCellShockPoints.at(ISHCELL);
                ISHCELLPOINT++) {
     m_X.at(ISHCELLPOINT) = X.at(ptrShCellPoint(ISHCELLPOINT,ISHCELL));
     m_Y.at(ISHCELLPOINT) = Y.at(ptrShCellPoint(ISHCELLPOINT,ISHCELL));
   }

   // interpolate the given coordinates
   // m_X and m_Y will be overwritten with the new values
   m_fitData(m_X,m_Y,SegPolynomialOrders.at(ISHCELL));

   nbCellShockPoints.at(ISHCELL)=m_X.size();

   // smooth the Y-coordinates
   if(SmoothingOption) { 
    cout << "        LeastSquaresMethod:: (!) smoothing option actived\n";
    iSmooth.curveForwardSmoothing(m_Y,6);
    m_Y = iSmooth.getYvector(); }

   // assign the new shock points coordinates to the vectors returned by 
   // the class
   for(unsigned ISHCELLPOINT=0;
                ISHCELLPOINT<nbCellShockPoints.at(ISHCELL);
                ISHCELLPOINT++) {
    Xnew.at(iShPointNew) = m_X.at(ISHCELLPOINT);
    Ynew.at(iShPointNew) = m_Y.at(ISHCELLPOINT);
    ++iShPointNew;
   }
  }

  Xnew.resize(iShPointNew);
  Ynew.resize(iShPointNew);

  // plot points interpolation
  FILE* fittedPoints;
  fittedPoints = fopen("log/LeastSquaresMethod.dat","w");

  fprintf(fittedPoints,"%s","TITLE = Fitted Points\n");
  fprintf(fittedPoints,"%s","VARIABLES = \"x0\" \"x1\"\n");
  fprintf(fittedPoints,"%s","ZONE T = \"Fitted points\"\n");
  fprintf(fittedPoints,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
  fprintf(fittedPoints,"%s %u %s","I= ",m_nbPoints,", J=1, K=1, ZONETYPE=Ordered\n");
  fprintf(fittedPoints,"%s","DATAPACKING=POINT\n");
  fprintf(fittedPoints,"%s","DT = (SINGLE, SINGLE)\n");

  for(unsigned IPOIN=0;IPOIN<m_nbPoints;IPOIN++) {
   fprintf(fittedPoints,"%22.14F",Xnew.at(IPOIN));
   fprintf(fittedPoints,"%22.14F",Ynew.at(IPOIN));
   fprintf(fittedPoints,"%s","\n");
  }

  fclose(fittedPoints);

  return;
}

//--------------------------------------------------------------------------//

