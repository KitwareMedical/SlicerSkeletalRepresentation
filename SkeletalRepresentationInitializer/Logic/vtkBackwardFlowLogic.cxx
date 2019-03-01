// This class provides logic of backward flow
// author: Zhiyuan Liu
// Date: Sept. 4, 2018
#include "vtkBackwardFlowLogic.h"
#include <iostream>
#include <string>

//#include "itkThinPlateSplineKernelTransform.h"
#include "itkThinPlateSplineExtended.h"
#include "itkPointSet.h"
//#include "itkTransformFileWriter.h"

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkMath.h"

#define EPS			1e-9
#define PI			3.141592653589793
#define EQZERO(x)	(fabs(x)<EPS)

void vtkBackwardFlowLogic::runApplyTPS()
{
    
}

void vtkBackwardFlowLogic::computePairwiseTPS(vtkPolyData* polyData_source, vtkPolyData* polyData_target, const char* outputFileName)
{
    typedef double CoordinateRepType;
//	typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType,3> TransformType;
	typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	// typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	// typedef PointSetType::Pointer PointSetPointer;
	typedef PointSetType::PointIdentifier PointIdType;

	PointSetType::Pointer sourceLandMarks = PointSetType::New();
	PointSetType::Pointer targetLandMarks = PointSetType::New();
	PointType p1; PointType p2; // same as double p1[3];
	PointSetType::PointsContainer::Pointer sourceLandMarkContainer
        = sourceLandMarks->GetPoints();
	PointSetType::PointsContainer::Pointer targetLandMarkContainer
        = targetLandMarks->GetPoints();

	PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
	PointIdType id_t = itk::NumericTraits< PointIdType >::Zero;
	// Read in the source points set
	for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i += 10){
		double p[3];
		polyData_source->GetPoint(i,p);
		p1[0] = p[0];
		p1[1] = p[1];
		p1[2] = p[2];
		sourceLandMarkContainer->InsertElement(id_s, p1);
		id_s++;
	}

	// Read in the target points set
	for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); i += 10){
		double p[3];
		polyData_target->GetPoint(i,p);
		p2[0] = p[0];
		p2[1] = p[1];
		p2[2] = p[2];
		targetLandMarkContainer->InsertElement(id_t, p2);
		id_t++;
	}

	TransformType::Pointer tps = TransformType::New();
	tps->SetSourceLandmarks(sourceLandMarks);
	tps->SetTargetLandmarks(targetLandMarks);

	cout<<"Computing W Matrix... "<<endl;
	tps->ComputeWMatrix();
	cout<<"Compute W Matrix finished!"<<endl;

//	PointType pos;
//	pos[0] = 0; pos[1] = 0; pos[2] = 0;
//	PointType new_pos = tps->TransformPoint(pos);

	// 1. write D matrix
	itkThinPlateSplineExtended::DMatrixType D = tps->getDMatrix();

	// 2. write A matrix
	itkThinPlateSplineExtended::AMatrixType A = tps->getAMatrix();

	// 3. write B vector
	itkThinPlateSplineExtended::BMatrixType B = tps->getBVector();

	std::ofstream fout;
	fout.open(outputFileName);

	if (fout)  {
//        for (unsigned int i=0;i< featurematrix.size();i++){
//            for (unsigned int j=0;j< featurematrix[i].size();j++)
//                fout<< featurematrix[i][j]<<" ";
//            fout<<endl;
//        }
		fout << D.rows() << "x" << D.cols() << "\n";
		for(unsigned int i = 0; i < D.rows(); i++) {
			for(unsigned int j = 0; j < D.cols(); j++) {
				fout << D.get(i,j)<<",";
			}
			fout << "\n";
		}

		fout << A.rows() << "x" << A.cols() << "\n";
		for(unsigned int i = 0; i < A.rows(); i++) {
			for(unsigned int j = 0; j < A.cols(); j++) {
				fout << A.get(i,j)<<",";
			}
			fout << "\n";
		}

		fout << B.size() << "\n";
		for(unsigned int i = 0; i < B.size(); i++) {
			fout << B.get(i)<<",";
		}
		fout << "\n";
		//cout<<"Successfully saved to: "<<filename<<endl;
	}
	else {
		cerr<<"Write out failed, cannot open the file!"<<endl;
		fout.flush();
		fout.close();
		return;
	}

	fout.close();

}
// void vtkBackwardFlowLogic::generateEllipsoidSrep(int numRow, int numCol, double ra, double rb, double rc, const char* outputPath)
// {
//    using namespace std;
//    ofstream fout(outputPath);
//    if (!fout.is_open()) {
//            cout << "File access failure!" <<endl;
//            return;
//    }

//    double mra = (ra*ra-rc*rc)/ra;
//    double mrb = (rb*rb-rc*rc)/rb;

//    fout << "model {" << endl
//             << "\tfigureCount = 1;" << endl
//             << "\tname = default;" << endl;
//    fout << "\tfigureTrees {" << endl
//             << "\t\tcount = 1;" << endl
//             << "\t\ttree[0] {" << endl
//             << "\t\t\tattachmentMode = 0;" << endl
//             << "\t\t\tblendAmount = 0;" << endl
//             << "\t\t\tblendExtent = 0;" << endl
//             << "\t\t\tchildCount = 0;" << endl
//             << "\t\t\tfigureId = 0;" << endl
//             << "\t\t\tlinkCount = 0;" << endl
//             << "\t\t}" << endl
//             << "\t}" << endl;
//    fout << "\tfigure[0] {" << endl
//             << "\t\tnumRows = " << numRow << ";" << endl
//             << "\t\tnumColumns = " << numCol << ";" << endl
//             << "\t\tnumLandmarks = 0;" << endl
//             << "\t\tpositivePolarity = 1;" << endl
//             << "\t\tpositiveSpace = 1;" << endl
//             << "\t\tsmoothness = 50;" << endl
//             << "\t\ttype = QuadFigure;" << endl;
//    fout << "\t\tcolor {" << endl
//             << "\t\t\tblue = 0;" << endl
//             << "\t\t\tgreen = 1;" << endl
//             << "\t\t\tred = 0;" << endl
//             << "\t\t}" << endl;

//    double ELLIPSE_SCALE = 0.9;
//    for (int row = 0; row < numRow; ++row) {
//            for (int col = 0; col < numCol; ++col) {
//                    fout << "\t\tprimitive[" << row << "][" << col << "] {" << endl;
//                    fout << "\t\t\tselected = 1;" << endl;
//                    if (row == 0 || row == numRow-1 || col == 0 || col == numCol-1) {
//                            fout << "\t\t\ttype = EndPrimitive;" << endl;
//                    } else {
//                            fout << "\t\t\ttype = StandardPrimitive;" << endl;
//                    }

//                    double x = 0;
//                    double y = 0;
//                    double z = 0;

//                    if (row == 1) {
//                            x = 2*mra*col / (numCol-1) - mra;
//                    } else {
//                            x = 2*mra*(col+1) / (numCol + 1) - mra;
//                    }
//                    y = (row - 1) * mrb*sqrt(1-x*x/(mra*mra));

//                    x*=ELLIPSE_SCALE;
//                    y*=ELLIPSE_SCALE;

////                    fout << "\t\t\tx = " << transform(x) << ";" << endl
////                             << "\t\t\ty = " << transform(y) << ";" << endl
////                             << "\t\t\tz = " << transform(z) << ";" << endl;

//                    double sinB = y*mra;
//                    double cosB = x*mrb;
//                    double l = normalize2(sinB, cosB);
//                    double cosA = l/(mra*mrb);
//                    double sinA = sqrt(1-cosA*cosA+EPS);

//                    double sx = ra * cosA * cosB - x;
//                    double sy = rb * cosA * sinB - y;
//                    double sz = rc * sinA ;

//                    double cx = ra * cosB - x;
//                    double cy = rb * sinB - y;
//                    double cz = 0;

//                    double len1 = normalize3(sx,sy,sz);
//                    double len2 = normalize2(cx,cy);
//                    if (EQZERO(len2)) cz = 1;

//                    fout << "\t\t\tux[0] = " << sx << ";" << endl
//                             << "\t\t\tux[1] = " << sx << ";" << endl
//                             << "\t\t\tux[2] = " << cx << ";" << endl
//                             << "\t\t\tuy[0] = " << sy << ";" << endl
//                             << "\t\t\tuy[1] = " << sy << ";" << endl
//                             << "\t\t\tuy[2] = " << cy << ";" << endl
//                             << "\t\t\tuz[0] = " <<-sz << ";" << endl
//                             << "\t\t\tuz[1] = " << sz << ";" << endl
//                             << "\t\t\tuz[2] = " << cz << ";" << endl;

//                    fout << "\t\t\tr[0] = " << scale(len1) << ";" << endl
//                             << "\t\t\tr[1] = " << scale(len1) << ";" << endl;
//                    if (row == 0 || row == numRow-1 || col == 0 || col == numCol-1) {
//                            fout << "\t\t\tr[2] = " << scale(len2) << ";" << endl;
//                    }

//                    fout << "\t\t}" << endl;
//            }
//    }

//    fout << "\t}" << endl;
//    fout << "\ttransformation {" << endl
//             << "\t\tscale = 1;" << endl
//             << "\t\trotation {" << endl
//             << "\t\t\tw = 1;" << endl
//             << "\t\t\tx = 0;" << endl
//             << "\t\t\ty = 0;" << endl
//             << "\t\t\tz = 0;" << endl
//             << "\t\t}" << endl
//             << "\t\ttranslation {" << endl
//             << "\t\t\tx = 0;" << endl
//             << "\t\t\ty = 0;" << endl
//             << "\t\t\tz = 0;" << endl
//             << "\t\t}" << endl
//             << "\t}" << endl;
//    fout << "}" << endl;
//    fout.flush();
//    fout.close();
// }
