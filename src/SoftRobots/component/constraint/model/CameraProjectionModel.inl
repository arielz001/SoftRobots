
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

#include <SoftRobots/component/constraint/model/CameraProjectionModel.h>
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
CameraProjectionModel<DataTypes>::CameraProjectionModel(MechanicalState* object)
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
    , d_radiusEllipse(initData(&d_radiusEllipse, "radiusEllipse",
                                 "Radius of the 3D feature ellipse/disk"))
                                 
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
CameraProjectionModel<DataTypes>::~CameraProjectionModel()
{
}





template<class DataTypes>
void CameraProjectionModel<DataTypes>::init()
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
void CameraProjectionModel<DataTypes>::reinit()
{
    internalInit();
}






template<class DataTypes>
void CameraProjectionModel<DataTypes>::internalInit()
{
    if(!d_radiusEllipse.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto radiusEllipse = sofa::helper::getReadAccessor(d_radiusEllipse);
    }
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


template<class DataTypes>
void CameraProjectionModel<DataTypes>::checkIndicesRegardingState()
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
void CameraProjectionModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int> > > defaultIndices = d_indices;
    defaultIndices.resize(1);
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}






Eigen::Matrix<double, 5, 1> calculateProjectedEllipse(
    double x, double y, double z,
    const Eigen::Matrix3d& R,
    double radius,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos)
{
    // relative coordinates of the camera 
    double x_rel = x - cameraPos[0];
    double y_rel = y - cameraPos[1];
    double z_rel = z - cameraPos[2];
    

    // 2d projection
    double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
    double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1]; 

    // inclination nad normal
    Eigen::Vector3d normal = R.col(2); 
    double cos_tilt = std::abs(normal(2));
    if (cos_tilt < 1e-3) cos_tilt = 1e-3;

    // semiaxes of the ellipse
    double semi_a = focalLength[0] * (radius / z_rel);  
    double semi_b = semi_a * cos_tilt;           

    // amgle of the semiaxes
    double alpha_rad = std::atan2(-normal(1), normal(0)) + (M_PI / 2.0);
    double alpha_deg = alpha_rad * (180.0 / M_PI);

    // normalization
    while (alpha_deg < 0.0) alpha_deg += 180.0;
    while (alpha_deg >= 180.0) alpha_deg -= 180.0;

    Eigen::Matrix<double, 5, 1> ellipse;
    ellipse << u, v, semi_a, semi_b, alpha_deg;

    std::cout << "Simulated Projection: " << ellipse.transpose() << std::endl;
    // std::cout << "relative position: " << x_rel << ", " << y_rel << ", " << z_rel << std::endl;
    // std::cout << "camera position: " << cameraPos[0] << ", " << cameraPos[1] << ", " << cameraPos[2] << std::endl;
    // std::cout << "xyz position: " << x << ", " << y << ", " << z << std::endl;

    return ellipse;
}





template<class DataTypes>
void CameraProjectionModel<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
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
void CameraProjectionModel<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
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
    
    // Accesores de SOFA para la cámara y el disco
    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();
    const double radius = d_radiusEllipse.getValue();
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); 


    auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
    const auto& readAccessor = sofa::helper::getReadAccessor(x);

    // epsilon in this case 
    const double Cambio = 1e-4;

    for (const auto& coord : readAccessor) {
        // get the 3d position 
        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // getting the orientation
        Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]); // (w, x, y, z)
        Eigen::Matrix3d R = q.toRotationMatrix();

        // this is the original projected ellipse without perturbation
        Eigen::Matrix<double, 5, 1> E_0 = calculateProjectedEllipse(x_pos, y_pos, z_pos, R, radius, focalLength, principalPoint, cameraPosition);

        // --- derivatives with respect to x ---
        Eigen::Matrix<double, 5, 1> E_x = calculateProjectedEllipse(x_pos + Cambio, y_pos, z_pos, R, radius, focalLength, principalPoint, cameraPosition);
        Eigen::Matrix<double, 5, 1> dE_dx = (E_x - E_0) / Cambio;

        // --- derivatives with respect to y ---
        Eigen::Matrix<double, 5, 1> E_y = calculateProjectedEllipse(x_pos, y_pos + Cambio, z_pos, R, radius, focalLength, principalPoint, cameraPosition);
        Eigen::Matrix<double, 5, 1> dE_dy = (E_y - E_0) / Cambio;

        // --- derivatives with respect to z ---
        Eigen::Matrix<double, 5, 1> E_z = calculateProjectedEllipse(x_pos, y_pos, z_pos + Cambio, R, radius, focalLength, principalPoint, cameraPosition);
        Eigen::Matrix<double, 5, 1> dE_dz = (E_z - E_0) / Cambio;

        // --- derivatives with respect to rotations (perturbation in X, Y, Z rotation) ---
        // perturbation Roll (around X)
        Eigen::Matrix3d R_rx = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitX());
        Eigen::Matrix<double, 5, 1> dE_dRx = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_rx, radius, focalLength, principalPoint, cameraPosition) - E_0) / Cambio;

        // perturbation Pitch (around Y)
        Eigen::Matrix3d R_ry = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitY());
        Eigen::Matrix<double, 5, 1> dE_dRy = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_ry, radius, focalLength, principalPoint, cameraPosition) - E_0) / Cambio;

        // perturbation Yaw (around Z)
        Eigen::Matrix3d R_rz = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitZ());
        Eigen::Matrix<double, 5, 1> dE_dRz = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_rz, radius, focalLength, principalPoint, cameraPosition) - E_0) / Cambio;

        // jacobian matrix (5 rows: u, v, a, b, alpha)
        for (int i = 0; i < 5; ++i) {
            Jacobian[i] = sofa::type::Vec<6, double>(
                dE_dx[i],   // d/dx
                dE_dy[i],   // d/dy
                dE_dz[i],   // d/dz
                dE_dRx[i],  // d/dRoll
                dE_dRy[i],  // d/dPitch
                dE_dRz[i]   // d/dYaw
            );
        } 


    }

    // write constraints in the global system of SOFA
    unsigned int index = 0;
    for (unsigned j = 0; j < 5; j++) { 
        MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
        rowIterator.setCol(0, Jacobian[j]);
        index++;
    }

    cIndex += index;
    cMatrix.endEdit();
    m_nbLines = cIndex - constraintIndex;
}




template<class DataTypes>


void CameraProjectionModel<DataTypes>::storeResults(vector<double> &delta)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    d_delta.setValue(delta);
}



template<class DataTypes>
void CameraProjectionModel<DataTypes>::setDefaultDirections()
{
    using Deriv = typename DataTypes::Deriv;
    sofa::type::vector<Deriv> jacobianInit(5); // 5 constraints (u, v, a, b, alpha)

    for (size_t i = 0; i < 5; ++i)
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
void CameraProjectionModel<DataTypes>::setDefaultUseDirections()
{
    Vec<Deriv::total_size, bool> useDirections;
    useDirections.assign(true);
    d_useDirections.setValue(useDirections);
}


template<class DataTypes>
void CameraProjectionModel<DataTypes>::normalizeDirections()
{

    WriteAccessor<sofa::Data<VecDeriv>> directions = d_directions;
    directions.resize(Deriv::total_size);
    for(unsigned int i=0; i<Deriv::total_size; i++)
        directions[i].normalize();
}



template<class DataTypes>
void CameraProjectionModel<DataTypes>::draw(const VisualParams* vparams)
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

    Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]);
    Eigen::Matrix3d R = q.toRotationMatrix();

    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();
    const double radius = d_radiusEllipse.getValue();
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); 
    // get the 5 parameters of the ellipse 2D [u, v, a, b, alpha]
    Eigen::Matrix<double, 5, 1> ellipse = calculateProjectedEllipse(
        x_pos, y_pos, z_pos, R, radius, focalLength, principalPoint, cameraPosition);

    double u = ellipse[0];
    double v = ellipse[1];
    double a = ellipse[2];
    double b = ellipse[3];
    double alpha = ellipse[4];

    // generate the contour of the ellipse (36 segments)
    const int num_segments = 36;
    sofa::type::vector<sofa::type::Vec3d> ellipsePoints3D;
    ellipsePoints3D.reserve(num_segments * 2);

    double cos_a = std::cos(alpha);
    double sin_a = std::sin(alpha);

    sofa::type::Vec3d prevPoint;

    for (int i = 0; i <= num_segments; ++i)
    {
        double theta = 2.0 * M_PI * i / num_segments;
        
        // parametric equation of the rotated ellipse in 2D (pixels)
        double x_local = a * std::cos(theta);
        double y_local = b * std::sin(theta);

        double u_p = u + (x_local * cos_a - y_local * sin_a);
        double v_p = v + (x_local * sin_a + y_local * cos_a);

        // projection back to 3D using the pinhole camera model
        double X_3d = (u_p - principalPoint[0]) * z_pos / focalLength[0];
        double Y_3d = (v_p - principalPoint[1]) * z_pos / focalLength[1];
        double Z_3d = z_pos;

        sofa::type::Vec3d currentPoint(X_3d, Y_3d, Z_3d);

        if (i > 0)
        {
            // add a point to form the line segment
            ellipsePoints3D.push_back(prevPoint);
            ellipsePoints3D.push_back(currentPoint);
        }
        prevPoint = currentPoint;
    }

    // render
    vparams->drawTool()->drawLines(ellipsePoints3D, 4.0f, RGBAColor::red());
    vector<Coord> points;
    points.push_back(positions[indices[0]]);
    drawPoints(vparams, points, 8.0f, RGBAColor::green());
}
} // namespace
