
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

using namespace Eigen;

namespace softrobots::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;
using sofa::core::ConstVecCoordId;
using sofa::helper::WriteAccessor;
using sofa::helper::ReadAccessor;
using sofa::type::vector;
using sofa::type::Vec;
using sofa::type::Vec3;
using sofa::type::RGBAColor;


template<class DataTypes>
CameraProjectionPointsModel<DataTypes>::CameraProjectionPointsModel(MechanicalState* object)
    : Inherit1(object)
    , d_indices(initData(&d_indices, "indices",
                                 "Indices of the point(s) to project"))

    , d_weight(initData(&d_weight, sofa::type::vector<Real>(Deriv::total_size, 1.), "weight",
                          "The parameter sets a weight to the minimization."))

    , d_directions(initData(&d_directions,"directions",
                          "Directions in which to solve position."))

    , d_Jacobian(initData(&d_Jacobian,"JacobianBq",
                          "Jacobian relating node motion to 2D image coordinates change"))

    , d_focalLength(initData(&d_focalLength, "focalLength",
                                 "Focal length of the camera (fx, fy)"))

    , d_principalPoint(initData(&d_principalPoint, "principalPoint",
                                 "Principal point of the camera (cx, cy)"))

    , d_useDirections(initData(&d_useDirections,"useDirections",
                              "Select directions to solve position."))

    , d_delta(initData(&d_delta, "delta", "Distance to target (2D residual)"))

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
                                    msg_info() << "Wrong size for weight data field. Resizing.";
                                    Real w = weight.empty() ? 1. : weight[0];
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

    if(m_state == nullptr)
    {
        msg_error() << "No mechanical state associated with this node. Deactivating object.";
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
    if(!d_directions.isSet())
        setDefaultDirections();
    else
        normalizeDirections();

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

    if(!d_indices.isSet())
    {
        msg_warning(this) << "Indices not defined. Default value assigned: 0.";
        setIndicesDefaultValue();
    }

    if(d_indices.getValue().size() > m_state->getSize())
    {
        msg_warning(this) << "Indices size exceeds points count. Resizing.";
        resizeIndicesRegardingState();
    }

    if(d_indices.getValue().size() == 0)
    {
        msg_error(this) << "Indices size is zero. Invalid component.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    checkIndicesRegardingState();
}


template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::checkIndicesRegardingState()
{
    ReadAccessor<sofa::Data<VecCoord>> positions = m_state->readPositions();

    if(d_indices.getValue().size() > positions.size())
    {
        msg_error(this) << "Indices size larger than mechanical state size";
        d_componentState = ComponentState::Invalid;
        return;
    }

    const auto& indices = d_indices.getValue();
    for(unsigned int i = 0; i < indices.size(); i++)
    {
        if (positions.size() <= indices[i])
        {
            msg_error(this) << "Index " << i << " exceeds mechanical state bounds";
            d_componentState = ComponentState::Invalid;
            return;
        }
    }
}


template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> defaultIndices = d_indices;
    defaultIndices.resize(1);
}

template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}


/**
 * @brief Proyecta un punto 3D al plano de imagen 2D (u, v) mediante el modelo Pinhole.
 */
Eigen::Vector2d calculateProjectedPoint(
    double x, double y, double z,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos)
{
    // Coordenadas relativas a la cámara
    double x_rel = x - cameraPos[0];
    double y_rel = y - cameraPos[1];
    double z_rel = z - cameraPos[2];

    // Evitar división por cero
    if (std::abs(z_rel) < 1e-6)
        z_rel = (z_rel >= 0) ? 1e-6 : -1e-6;

    // Proyección 2D Pinhole
    double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
    double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1];

    return Eigen::Vector2d(u, v);
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

    const double Cambio = 1e-4; // Paso de perturbación

    for (const auto& coord : readAccessor) {
        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // Proyección original en 2D (u, v)
        Eigen::Vector2d P_0 = calculateProjectedPoint(x_pos, y_pos, z_pos, focalLength, principalPoint, cameraPosition);

        // Derivadas numéricas respecto a x, y, z
        Eigen::Vector2d P_x = calculateProjectedPoint(x_pos + Cambio, y_pos, z_pos, focalLength, principalPoint, cameraPosition);
        Eigen::Vector2d dP_dx = (P_x - P_0) / Cambio;

        Eigen::Vector2d P_y = calculateProjectedPoint(x_pos, y_pos + Cambio, z_pos, focalLength, principalPoint, cameraPosition);
        Eigen::Vector2d dP_dy = (P_y - P_0) / Cambio;

        Eigen::Vector2d P_z = calculateProjectedPoint(x_pos, y_pos, z_pos + Cambio, focalLength, principalPoint, cameraPosition);
        Eigen::Vector2d dP_dz = (P_z - P_0) / Cambio;

        // Matriz Jacobiana de dimensión 2 x Deriv::total_size (u y v)
        for (int i = 0; i < 2; ++i) {
            // Nota: Se asume que Deriv maneja al menos 3 Grados de Libertad traslacionales
            sofa::type::Vec<6, double> row(0.0);
            row[0] = dP_dx[i];
            row[1] = dP_dy[i];
            row[2] = dP_dz[i];
            // Componentes de rotación quedan en 0 si las hay
            Jacobian[i] = row;
        } 
    }

    // Escribir restricciones globales (2 restricciones: u, v)
    unsigned int index = 0;
    for (unsigned j = 0; j < 2; j++) { 
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
    sofa::type::vector<Deriv> jacobianInit(2); // 2 restricciones (u, v)

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
    for(unsigned int i = 0; i < Deriv::total_size; i++)
        directions[i].normalize();
}


template<class DataTypes>
void CameraProjectionPointsModel<DataTypes>::draw(const VisualParams* vparams)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    ReadAccessor<sofa::Data<VecCoord>> positions = m_state->readPositions();
    ReadAccessor<sofa::Data<sofa::type::vector<sofa::Index>>> indices = d_indices;

    if (indices.empty() || positions.empty()) 
        return;

    // Posición 3D del punto seleccionado
    const auto& coord = positions[indices[0]];
    double x_pos = coord[0];
    double y_pos = coord[1];
    double z_pos = coord[2];

    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); 

    // Obtener proyección 2D (u, v)
    Eigen::Vector2d point2D = calculateProjectedPoint(
        x_pos, y_pos, z_pos, focalLength, principalPoint, cameraPosition);

    // Re-proyectar al espacio 3D a la profundidad z_pos para visualización en escena SOFA
    double X_3d = (point2D[0] - principalPoint[0]) * z_pos / focalLength[0];
    double Y_3d = (point2D[1] - principalPoint[1]) * z_pos / focalLength[1];
    double Z_3d = z_pos;

    vector<Coord> points;
    points.push_back(Coord(X_3d, Y_3d, Z_3d));

    // Dibujar el punto proyectado en color rojo
    drawPoints(vparams, points, 8.0f, RGBAColor::red());
}

} // namespace softrobots::constraint