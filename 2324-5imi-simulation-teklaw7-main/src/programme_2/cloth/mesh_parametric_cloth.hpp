/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifndef MESH_PARAMETRIC_CLOTH_HPP
#define MESH_PARAMETRIC_CLOTH_HPP

#include "../lib/mesh/mesh_parametric.hpp"
#include "../lib/common/exception_cpe.hpp"

namespace cpe
{
class mesh_parametric_cloth : public mesh_parametric
{
public:
    using mesh_parametric::mesh_parametric;

    void set_plane_xy_unit(int const size_u_param,int const size_v_param);

    vec3 const& speed(int ku,int kv) const;
    vec3& speed(int ku,int kv);

    vec3 const& force(int ku,int kv) const;
    vec3& force(int ku,int kv);

    void update_force();
    void integration_step(float dt);
    vec3 spring_compute_force(int u1, int v1, int u2, int v2, float const L0, float const K);
    void intersect_plane(int ku, int kv);
    void intersect_sphere(int ku, int kv, vec3 const& center, float const radius);

private:

    std::vector<vec3> speed_data;
    std::vector<vec3> force_data;

};

class exception_divergence : public exception_cpe
{
    using exception_cpe::exception_cpe;
};

}

#endif
