/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
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

/*! \file Case_RayleighPlessetBubble.hpp
    \brief Header for Case_RayleighPlessetBubble.cpp
*/

/*! \class Case_RayleighPlessetBubble
    \brief Derived class: handles the dynamics of bubbles that have been initiated based on the Rayleigh-Plesset model.

    The bubble initial condition are taken from the Rayleigh-Plesset model with R0 computed based on
    epsilon and k, and the potentials on the bubble taken as 0.
    This derived class inherits members from BubbleData.cpp and BubbleData.hpp.
*/

#ifndef CASE_RAYLEIGHPLESSETBUBBLE_HPP
#define CASE_RAYLEIGHPLESSETBUBBLE_HPP

#include "BubbleData.hpp"

class Case_RayleighPlessetBubble : public BubbleData {
public:
    template<typename Input>
    Case_RayleighPlessetBubble(Input &data): BubbleData(data) {}

    virtual ~Case_RayleighPlessetBubble() {}

    /*! \brief Initializes the node points on the bubble surface.*/
    void initialize() override;

};


#endif // CASE_RAYLEIGHPLESSETBUBBLE_HPP
