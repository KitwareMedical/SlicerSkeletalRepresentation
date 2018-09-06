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

void vtkBackwardFlowLogic::runApplyTPS()
{
    
}

void vtkBackwardFlowLogic::computePairwiseTPS(vtkPolyData* polyData_source, vtkPolyData* polyData_target, const char* outputFileName)
{
    typedef double CoordinateRepType;
//	typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType,3> TransformType;
	typedef itkThinPlateSplineExtended TransformType;
	typedef itk::Point< CoordinateRepType, 3 > PointType;
	typedef std::vector< PointType > PointArrayType;
	typedef TransformType::PointSetType PointSetType;
	typedef PointSetType::Pointer PointSetPointer;
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
//        for (int i=0;i< featurematrix.size();i++){
//            for (int j=0;j< featurematrix[i].size();j++)
//                fout<< featurematrix[i][j]<<" ";
//            fout<<endl;
//        }
		fout << D.rows() << "x" << D.cols() << "\n";
		for(int i = 0; i < D.rows(); i++) {
			for(int j = 0; j < D.cols(); j++) {
				fout << D.get(i,j)<<",";
			}
			fout << "\n";
		}

		fout << A.rows() << "x" << A.cols() << "\n";
		for(int i = 0; i < A.rows(); i++) {
			for(int j = 0; j < A.cols(); j++) {
				fout << A.get(i,j)<<",";
			}
			fout << "\n";
		}

		fout << B.size() << "\n";
		for(int i = 0; i < B.size(); i++) {
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