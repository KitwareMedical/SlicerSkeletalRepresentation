/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/
//
// Created by jphong on 2/19/18.
//

#ifndef __itkThinPlateSplineExtended_H__
#define __itkThinPlateSplineExtended_H__

#include "itkThinPlateSplineKernelTransform.h"


class itkThinPlateSplineExtended :
		public itk::ThinPlateSplineKernelTransform<double, 3> {
public:
	/** Standard class typedefs. */
	typedef itkThinPlateSplineExtended                     Self;
	typedef itk::ThinPlateSplineKernelTransform<double, 3> Superclass;
	typedef itk::SmartPointer<Self>                            Pointer;
	typedef itk::SmartPointer<const Self>                      ConstPointer;

	/** New macro for creation of through a Smart Pointer */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( itkThinPlateSplineExtended, KernelTransform );


//	/** The Deformation matrix.
//		This is an auxiliary matrix that will hold the
//		Deformation (non-affine) part of the transform.
//		Those are the coefficients that will multiply the
//		Kernel function */
//	DMatrixType m_DMatrix;
	DMatrixType getDMatrix() {return m_DMatrix;};
	void setDMatrix(DMatrixType D) {m_DMatrix = D;};
//
//	/** Rotatinoal/Shearing part of the Affine component of the Transformation */
//	AMatrixType m_AMatrix;
	AMatrixType getAMatrix() {return m_AMatrix;};
	void setAMatrix(AMatrixType A) {m_AMatrix = A;};
//
//	/** Translational part of the Affine component of the Transformation */
//	BMatrixType m_BVector;
	BMatrixType getBVector() {return m_BVector;};
	void setBVector(BMatrixType B) {m_BVector = B;};
protected:
	itkThinPlateSplineExtended() {};
	virtual ~itkThinPlateSplineExtended() {}
private:
	itkThinPlateSplineExtended(const Self&); //purposely not implemented
	void operator=(const Self&); //purposely not implemented
};


#endif //__itkThinPlateSplineExtended_H__
