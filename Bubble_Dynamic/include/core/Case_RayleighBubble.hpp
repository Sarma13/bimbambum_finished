/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * Case_RayleighBubble.hpp
 *
 * Description:
 * Header of the derived class containing function definitions, properties and
 * quantities associated with the bubble dynamics. Here, the bubble initial
 * condition are taken from the Rayleigh model with R0 = 0.1 and the potentials
 * on the bubble surface derived from the Rayleigh equation.
 * This class inherits from BubbleData.cpp/BubbleData.hpp.
 *
 * -----------------------------------------------------------------------------
 * Copyright (C) 2023 Armand Sieber
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file Case_RayleighBubble.hpp
    \brief Header for Case_RayleighBubble.cpp
*/

/*! \class Case_RayleighBubble
    \brief Derived class: handles the dynamics of bubbles that have been initiated based on the Rayleigh model.

    The bubble initial condition are taken from the Rayleigh model with R0 = 0.1 and the potentials
    on the bubble surface derived from the Rayleigh equation.
    This derived class inherits members from BubbleData.cpp and BubbleData.hpp.
*/

#ifndef CPP_CODE_CASE_RAYLEIGHBUBBLE_HPP
#define CPP_CODE_CASE_RAYLEIGHBUBBLE_HPP

#include "BubbleData.hpp"

class Case_RayleighBubble : public BubbleData {
public:
    template<typename Input>
    Case_RayleighBubble(Input &data): BubbleData(data) {}

    virtual ~Case_RayleighBubble() {}

    /*! \brief Initializes the node points on the bubble surface.*/
    void initialize() override;

};


#endif //CPP_CODE_CASE_RAYLEIGHBUBBLE_H
