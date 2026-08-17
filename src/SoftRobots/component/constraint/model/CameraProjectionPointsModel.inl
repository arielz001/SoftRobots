
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
*                           Plugin SoftRobots                                 *
*                                                                             *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
******************************************************************************/
#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <SoftRobots/component/constraint/model/CameraProjectionPointsModel.h>
#include <Eigen/Dense>

//#include <Eigen/Geometry>  // Para cuaterniones
//using namespace std;
using namespace Eigen;

namespace softrobots::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;
using sofa::core::ConstVecCoordId ;
using sofa::helper::WriteAccessor ;
using sofa::helper::ReadAccessor ;
using sofa::type::vector ;
using sofa::type::Vec;
using sofa::type::Vec3;
using sofa::type::RGBAColor;



template<class DataTypes>
CameraProjectionPointsModel<DataTypes>::CameraProjectionPointsModel(MechanicalState* object)
    : Inherit1(object)
    , d_indices(initData(&d_indices, "indices",
                                 "If indices size is lower than target size, \n"
                                 "some target will not be considered"))

    , d_weight(initData(&d_weight, sofa::type::vector<Real>(Deriv::total_size, 1.), "weight",
                          "The parameter sets a weight to the minimization."))

    , d_directions(initData(&d_directions,"directions",
                          "The parameter directions allows to specify the directions in \n"
                          "which you want to solve the position."))

    , d_Jacobian(initData(&d_Jacobian,"JacobianBq",
                          "Jacobian relating node motion to B change\n"
                          "."))
                                 
    , d_focalLength(initData(&d_focalLength, "focalLength",
                                 "Focal length of the camera"))

    , d_principalPoint(initData(&d_principalPoint, "principalPoint",
                                 "Principal point of the camera"))

    , d_useDirections(initData(&d_useDirections,"useDirections",
                              "The parameter useDirections allows to select the directions in \n"
                              "which you want to solve the position. If unspecified, the default \n"
                              "values are all true."))

    , d_delta(initData(&d_delta, "delta","Distance to target"))

    , d_cameraPosition(initData(&d_cameraPosition, "cameraPosition",
                                 "Position of the camera in the scene"))
{
    d_delta.setReadOnly(true);

    this->addUpdateCallback("updateWeight", {&d_weight}, [this](const sofa::core::DataTracker& t)
                            {
                                SOFA_UNUSED(t);
                                auto weight = sofa::helper::getWriteAccessor(d_weight);
                                if (weight.size() != Deriv::total_size)
                                {
                                    msg_info() << "Wrong size for the data field weight, " << weight.size() <<
                                        " instead of " << Deriv::total_size << ". Resizing, with weight[0] as the default value.";
                                    Real w = weight.empty()? 1.: weight[0];
                                    d_weight.setValue(sofa::type::vector<Real>(Deriv::total_size, w));
                                }
                                return sofa::core::objectmodel::ComponentState::Valid;
                            }, {});
}




template<class DataTypes>
CameraProjectionPointsModel<DataTypes>::~CameraProjectionPointsModel()
{
}





template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::init()
{
    d_componentState = ComponentState::Valid;
    Inherit1::init();

    if(m_state==nullptr)
    {
        msg_error() << "There is no mechanical state associated with this node. "
                        "the object is deactivated. "
                        "To remove this error message fix your scene possibly by "
                        "adding a MechanicalObject." ;
        d_componentState = ComponentState::Invalid;
        return;
    }

    internalInit();
}





template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::reinit()
{
    internalInit();
}






template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::internalInit()
{
    if(!d_focalLength.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        // const auto mum = sofa::helper::getReadAccessor(d_mum);
        const auto focalLength = sofa::helper::getReadAccessor(d_focalLength);
    }
    if(!d_principalPoint.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto principalPoint = sofa::helper::getReadAccessor(d_principalPoint);
    }

    if(!d_cameraPosition.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto cameraPosition = sofa::helper::getReadAccessor(d_cameraPosition);
    }


// ############################
    if(!d_directions.isSet())
        setDefaultDirections();
    else
        normalizeDirections();
// ###########################


    if(!d_useDirections.isSet())
    {
        setDefaultUseDirections();
    }
    else
    {
        const auto useDirections = sofa::helper::getReadAccessor(d_useDirections);
        if (std::find(useDirections.begin(), useDirections.end(), true) == useDirections.end())
        {
            setDefaultUseDirections();
            msg_warning(this) << "No direction given in useDirection. Set default all.";
        }
    }
// #########################
    if(!d_indices.isSet())
    {
        msg_warning(this) <<"Indices not defined. Default value assigned 0.";
        setIndicesDefaultValue();
    }

    if(d_indices.getValue().size() > m_state->getSize())
    {
        msg_warning(this) <<"Indices size can not be larger than the number of point in the context. Launch resize process.";
        resizeIndicesRegardingState();
    }

    if(d_indices.getValue().size() == 0)
    {
        msg_error(this) <<"Indices size is zero. The component will not work.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    checkIndicesRegardingState();
}



// =============================================================================

template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::checkIndicesRegardingState()
{
    ReadAccessor<sofa::Data<VecCoord> > positions = m_state->readPositions();

    if(d_indices.getValue().size() > positions.size())
    {
        msg_error(this) << "Indices size is larger than mechanicalState size" ;
        d_componentState = ComponentState::Invalid;
        return;
    }

    const auto& indices = d_indices.getValue();
    for(unsigned int i=0; i<indices .size(); i++)
    {
        if (positions.size() <= indices[i])
        {
            msg_error(this) << "Index at index " << i << " is too large regarding mechanicalState [position] size" ;
            d_componentState = ComponentState::Invalid;
            return;
        }
    }
}




template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int> > > defaultIndices = d_indices;
    defaultIndices.resize(1);
}

template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}





Eigen::Matrix<double, 2, 1> calculateProjectedPoint(
    double x, double y, double z,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos)
{
    // relative coordinates of the camera 
    double x_rel = x - cameraPos[0];
    double y_rel = y - cameraPos[1];
    double z_rel = z - cameraPos[2];
    // std::cout << "z_rel: " << z_rel << std::endl;
    // 2d projection
    double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
    double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1]; 


    Eigen::Matrix<double, 2, 1> point;
    point << u, v;

    return point;
}



template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                               sofa::linearalgebra::BaseVector *resV,
                                                               const sofa::linearalgebra::BaseVector *Jdx)
{
    SOFA_UNUSED(cParams);
    SOFA_UNUSED(Jdx);

    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
    const auto& delta = sofa::helper::getReadAccessor(d_delta);

    for (size_t i = 0; i < delta.size(); ++i)
    {
        resV->set(constraintIndex + i, delta[i]);
    }
}




template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
                                                              DataMatrixDeriv &cMatrix,
                                                              unsigned int &cIndex,
                                                              const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    d_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
    MatrixDeriv& column = *cMatrix.beginEdit();
    
    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();

    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); 

    auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
    const auto& readAccessor = sofa::helper::getReadAccessor(x);

    // epsilon in this case 
    const double Cambio = 1e-4;

    // for (const auto& coord : readAccessor) {
    //     // get the 3d position 
    //     double x_pos = coord[0];
    //     double y_pos = coord[1];
    //     double z_pos = coord[2];


    //     // this is the original projected point without perturbation
    //     Eigen::Matrix<double, 2, 1> E_0 = calculateProjectedPoint(x_pos, y_pos, z_pos, focalLength, principalPoint, cameraPosition);

    //     // --- derivatives with respect to x ---
    //     Eigen::Matrix<double, 2, 1> E_x = calculateProjectedPoint(x_pos + Cambio, y_pos, z_pos, focalLength, principalPoint, cameraPosition);
    //     Eigen::Matrix<double, 2, 1> dE_dx = (E_x - E_0) / Cambio;

    //     // --- derivatives with respect to y ---
    //     Eigen::Matrix<double, 2, 1> E_y = calculateProjectedPoint(x_pos, y_pos + Cambio, z_pos, focalLength, principalPoint, cameraPosition);
    //     Eigen::Matrix<double, 2, 1> dE_dy = (E_y - E_0) / Cambio;

    //     // --- derivatives with respect to z ---
    //     Eigen::Matrix<double, 2, 1> E_z = calculateProjectedPoint(x_pos, y_pos, z_pos + Cambio, focalLength, principalPoint, cameraPosition);
    //     Eigen::Matrix<double, 2, 1> dE_dz = (E_z - E_0) / Cambio;


    //     // jacobian matrix (2 rows: u, v)
    //     for (int i = 0; i < 2; ++i) {
    //         Jacobian[i] = sofa::type::Vec<3, double>(
    //             dE_dx[i],   // d/dx
    //             dE_dy[i],   // d/dy
    //             // 0.0   // d/dz
    //             dE_dz[i]   // d/dz
    //         );
    //     } 


    // }

    for (size_t pointIdx = 0; pointIdx < readAccessor.size(); ++pointIdx) {
        const auto& coord = readAccessor[pointIdx];

        const double X = coord[0] - cameraPosition[0];
        const double Y = coord[1] - cameraPosition[1];
        double Z = coord[2] - cameraPosition[2];

        if (std::abs(Z) < 1e-6) {
            Z = (Z >= 0) ? 1e-6 : -1e-6;
        }

        const double invZ = 1.0 / Z;
        const double invZ2 = invZ * invZ;

        
        const double fx = focalLength[0];
        const double fy = (focalLength.size() > 1) ? focalLength[1] : focalLength[0];

        // [du/dx, du/dy, du/dz]
        sofa::type::Vec<3, double> dU( fx * invZ, 0.0, -fx * X * invZ2 );

        // [dv/dx, dv/dy, dv/dz]
        sofa::type::Vec<3, double> dV(0.0,  fy * invZ, -fy * Y * invZ2 );

        Jacobian[2 * pointIdx]     = dU;
        Jacobian[2 * pointIdx + 1] = dV;
    }

    // write constraints in the global system of SOFA
    unsigned int index = 0;
    for (unsigned j = 0; j < 3; j++) { 
        MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
        rowIterator.setCol(0, Jacobian[j]);
        index++;
    }

    cIndex += index;
    cMatrix.endEdit();
    m_nbLines = cIndex - constraintIndex;
}




template<class DataTypes>


void CameraProjectionPointsModel<DataTypes>::storeResults(vector<double> &delta)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    d_delta.setValue(delta);
}



template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::setDefaultDirections()
{
    using Deriv = typename DataTypes::Deriv;
    sofa::type::vector<Deriv> jacobianInit(2); // 2 constraints (u, v)

    for (size_t i = 0; i < 2; ++i)
    {
        jacobianInit[i].clear();
        if (i < Deriv::total_size)
        {
            jacobianInit[i][i] = 1.0;
        }
    }
    d_Jacobian.setValue(jacobianInit);
}

template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::setDefaultUseDirections()
{
    Vec<Deriv::total_size, bool> useDirections;
    useDirections.assign(true);
    d_useDirections.setValue(useDirections);
}


template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::normalizeDirections()
{

    WriteAccessor<sofa::Data<VecDeriv>> directions = d_directions;
    directions.resize(Deriv::total_size);
    for(unsigned int i=0; i<Deriv::total_size; i++)
        directions[i].normalize();
}


template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::draw(const VisualParams* vparams)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    // this is for drawing only if the option to show constraints/forcefields is activated
    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    ReadAccessor<sofa::Data<VecCoord>> positions = m_state->readPositions();
    ReadAccessor<sofa::Data<sofa::type::vector<sofa::Index>>> indices = d_indices;

    if (indices.empty() || positions.empty()) 
        return;

    // get 3d position
    const auto& coord = positions[indices[0]];
    double x_pos = coord[0];
    double y_pos = coord[1];
    double z_pos = coord[2];


    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();

    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); 
    // get the 2 parameters of the point 2D [u, v]
    Eigen::Matrix<double, 2, 1> point = calculateProjectedPoint(
        x_pos, y_pos, z_pos, focalLength, principalPoint, cameraPosition);

    double u = point[0];
    double v = point[1];

    double z_rel = z_pos - cameraPosition[2];

    double X_3d = (u - principalPoint[0]) * z_rel / focalLength[0] + cameraPosition[0];
    double Y_3d = (v - principalPoint[1]) * z_rel / focalLength[1] + cameraPosition[1];
    double Z_3d = z_pos;
    // --------------------------------------------------


    sofa::type::Vec3d projectedPoint3D(X_3d, Y_3d, Z_3d);

    sofa::type::vector<sofa::type::Vec3d> linePoints;
    linePoints.push_back(cameraPosition);     
    linePoints.push_back(projectedPoint3D);   

    vparams->drawTool()->drawLines(linePoints, 4.0f,  RGBAColor::red());

    vector<Coord> points;
    points.push_back(positions[indices[0]]);
    drawPoints(vparams, points, 8.0f, RGBAColor::green());
}
}


