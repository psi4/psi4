/*-
 * Copyright (c) 2012-2014 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include <math_util.h>

#include "common.h"
#include "rand.h"

#define MAX_ITER 10

struct body {
	mat_t rotmat;
	vec_t pos;
	vec_t vel;
	vec_t vel_old;
	vec_t angmom;
	vec_t angmom_old;
	vec_t force;
	vec_t torque;
	vec_t inertia;
	vec_t inertia_inv;
	double mass;
};

struct nvt_data {
	double chi;
	double chi_dt;
};

struct npt_data {
	double chi;
	double chi_dt;
	double eta;
};

struct md {
	size_t n_bodies;
	struct body *bodies;
	size_t n_freedom;
	vec_t box;
	double potential_energy;
	double (*get_invariant)(const struct md *);
	void (*update_step)(struct md *);
	struct state *state;
	void *data;
};

void sim_md(struct state *state);

static vec_t wrap(const struct md *md, const vec_t *pos)
{
	if (!cfg_get_bool(md->state->cfg, "enable_pbc"))
		return *pos;

	vec_t sub = {
		md->box.x * floor(pos->x / md->box.x),
		md->box.y * floor(pos->y / md->box.y),
		md->box.z * floor(pos->z / md->box.z)
	};

	return vec_sub(pos, &sub);
}

static double get_kinetic_energy(const struct md *md)
{
	double ke = 0.0;

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		ke += body->mass * body->vel.x * body->vel.x;
		ke += body->mass * body->vel.y * body->vel.y;
		ke += body->mass * body->vel.z * body->vel.z;

		ke += body->angmom.x * body->angmom.x * body->inertia_inv.x;
		ke += body->angmom.y * body->angmom.y * body->inertia_inv.y;
		ke += body->angmom.z * body->angmom.z * body->inertia_inv.z;
	}

	return 0.5 * ke;
}

static double get_temperature(const struct md *md)
{
	double ke = get_kinetic_energy(md);

	return 2.0 * ke / BOLTZMANN / md->n_freedom;
}

static double get_volume(const struct md *md)
{
	return md->box.x * md->box.y * md->box.z;
}

static double get_pressure(const struct md *md)
{
	double volume = get_volume(md);
	vec_t pressure = vec_zero;

	for (size_t i = 0; i < md->n_bodies; i++) {
		const struct body *body = md->bodies + i;

		pressure.x += body->mass * body->vel.x * body->vel.x;
		pressure.y += body->mass * body->vel.y * body->vel.y;
		pressure.z += body->mass * body->vel.z * body->vel.z;
	}

	mat_t stress;
	check_fail(efp_get_stress_tensor(md->state->efp, (double *)&stress));

	pressure.x = (pressure.x + stress.xx) / volume;
	pressure.y = (pressure.y + stress.yy) / volume;
	pressure.z = (pressure.z + stress.zz) / volume;

	return (pressure.x + pressure.y + pressure.z) / 3.0;
}

static double get_invariant_nve(const struct md *md)
{
	return md->potential_energy + get_kinetic_energy(md);
}

static double get_invariant_nvt(const struct md *md)
{
	struct nvt_data *data = (struct nvt_data *)md->data;

	double t_tau = cfg_get_double(md->state->cfg, "thermostat_tau");
	double t_target = cfg_get_double(md->state->cfg, "temperature");

	double kbt = BOLTZMANN * t_target;

	double t_virt = kbt * md->n_freedom * (data->chi_dt +
				data->chi * data->chi * t_tau * t_tau / 2.0);

	return md->potential_energy + get_kinetic_energy(md) + t_virt;
}

static double get_invariant_npt(const struct md *md)
{
	struct npt_data *data = (struct npt_data *)md->data;

	double t_tau = cfg_get_double(md->state->cfg, "thermostat_tau");
	double t_target = cfg_get_double(md->state->cfg, "temperature");
	double p_tau = cfg_get_double(md->state->cfg, "barostat_tau");
	double p_target = cfg_get_double(md->state->cfg, "pressure");

	double kbt = BOLTZMANN * t_target;
	double volume = get_volume(md);

	double t_virt = kbt * md->n_freedom * (data->chi_dt +
				data->chi * data->chi * t_tau * t_tau / 2.0);

	double p_virt = p_target * volume + 3.0 * md->n_bodies * kbt *
				data->eta * data->eta * p_tau * p_tau / 2.0;

	return md->potential_energy + get_kinetic_energy(md) + t_virt + p_virt;
}

static vec_t get_system_com(const struct md *md)
{
	double mass = 0.0;
	vec_t com = vec_zero;

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;
		vec_t pos = wrap(md, &body->pos);

		com.x += body->mass * pos.x;
		com.y += body->mass * pos.y;
		com.z += body->mass * pos.z;

		mass += body->mass;
	}

	com.x /= mass;
	com.y /= mass;
	com.z /= mass;

	return com;
}

static vec_t get_system_com_velocity(const struct md *md)
{
	double mass = 0.0;
	vec_t cv = vec_zero;

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		cv.x += body->vel.x * body->mass;
		cv.y += body->vel.y * body->mass;
		cv.z += body->vel.z * body->mass;

		mass += body->mass;
	}

	cv.x /= mass;
	cv.y /= mass;
	cv.z /= mass;

	return cv;
}

static vec_t get_system_angular_momentum(const struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);

	vec_t am = vec_zero;

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t pos = wrap(md, &body->pos);
		vec_t dr = vec_sub(&pos, &cp);
		vec_t dv = vec_sub(&body->vel, &cv);

		am.x += (dr.y * dv.z - dr.z * dv.y) * body->mass;
		am.y += (dr.z * dv.x - dr.x * dv.z) * body->mass;
		am.z += (dr.x * dv.y - dr.y * dv.x) * body->mass;
	}

	return am;
}

static mat_t get_system_inertia_tensor(const struct md *md)
{
	mat_t inertia = mat_zero;
	vec_t com = get_system_com(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t pos = wrap(md, &body->pos);
		vec_t dr = vec_sub(&pos, &com);

		inertia.xx += body->mass * (dr.y * dr.y + dr.z * dr.z);
		inertia.yy += body->mass * (dr.x * dr.x + dr.z * dr.z);
		inertia.zz += body->mass * (dr.x * dr.x + dr.y * dr.y);
		inertia.xy -= body->mass * dr.x * dr.y;
		inertia.xz -= body->mass * dr.x * dr.z;
		inertia.yz -= body->mass * dr.y * dr.z;
	}

	inertia.yx = inertia.xy;
	inertia.zx = inertia.xz;
	inertia.zy = inertia.yz;

	return inertia;
}

static void remove_system_drift(struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);
	vec_t am = get_system_angular_momentum(md);

	mat_t inertia = get_system_inertia_tensor(md);
	mat_t inertia_inv = mat_zero;

	if (inertia.xx < EPSILON ||
	    inertia.yy < EPSILON ||
	    inertia.zz < EPSILON) {
		inertia_inv.xx = inertia.xx < EPSILON ? 0.0 : 1.0 / inertia.xx;
		inertia_inv.yy = inertia.yy < EPSILON ? 0.0 : 1.0 / inertia.yy;
		inertia_inv.zz = inertia.zz < EPSILON ? 0.0 : 1.0 / inertia.zz;
	}
	else {
		inertia_inv = mat_inv(&inertia);
	}

	vec_t av = mat_vec(&inertia_inv, &am);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;
		vec_t pos = wrap(md, &body->pos);

		vec_t cross = {
			av.y * (pos.z - cp.z) - av.z * (pos.y - cp.y),
			av.z * (pos.x - cp.x) - av.x * (pos.z - cp.z),
			av.x * (pos.y - cp.y) - av.y * (pos.x - cp.x)
		};

		body->vel.x -= cv.x + cross.x;
		body->vel.y -= cv.y + cross.y;
		body->vel.z -= cv.z + cross.z;
	}

	vec_t cv2 = get_system_com_velocity(md);
	vec_t am2 = get_system_angular_momentum(md);

	assert(vec_len(&cv2) < EPSILON && vec_len(&am2) < EPSILON);
}

static void compute_forces(struct md *md)
{
	for (size_t i = 0; i < md->n_bodies; i++) {
		double crd[12];

		memcpy(crd, &md->bodies[i].pos, 3 * sizeof(double));
		memcpy(crd + 3, &md->bodies[i].rotmat, 9 * sizeof(double));

		check_fail(efp_set_frag_coordinates(md->state->efp, i, EFP_COORD_TYPE_ROTMAT, crd));
	}

	compute_energy(md->state, true);

	md->potential_energy = md->state->energy;

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->force.x = -md->state->grad[6 * i + 0];
		body->force.y = -md->state->grad[6 * i + 1];
		body->force.z = -md->state->grad[6 * i + 2];

		body->torque.x = -md->state->grad[6 * i + 3];
		body->torque.y = -md->state->grad[6 * i + 4];
		body->torque.z = -md->state->grad[6 * i + 5];

		/* convert torque to body frame */
		body->torque = mat_trans_vec(&body->rotmat, &body->torque);
	}
}

static void set_body_mass_and_inertia(struct efp *efp, size_t idx, struct body *body)
{
	double mass, inertia[3];

	check_fail(efp_get_frag_mass(efp, idx, &mass));
	check_fail(efp_get_frag_inertia(efp, idx, inertia));

	body->mass = AMU_TO_AU * mass;

	body->inertia.x = AMU_TO_AU * inertia[0];
	body->inertia.y = AMU_TO_AU * inertia[1];
	body->inertia.z = AMU_TO_AU * inertia[2];

	body->inertia_inv.x = body->inertia.x < EPSILON ? 0.0 : 1.0 / body->inertia.x;
	body->inertia_inv.y = body->inertia.y < EPSILON ? 0.0 : 1.0 / body->inertia.y;
	body->inertia_inv.z = body->inertia.z < EPSILON ? 0.0 : 1.0 / body->inertia.z;
}

static void rotate_step(size_t a1, size_t a2, double angle, vec_t *angmom, mat_t *rotmat)
{
	mat_t rot = { 1.0, 0.0, 0.0,
		      0.0, 1.0, 0.0,
		      0.0, 0.0, 1.0 };

	double cosa = cos(angle);
	double sina = sin(angle);

	mat_set(&rot, a1, a1,  cosa);
	mat_set(&rot, a2, a2,  cosa);
	mat_set(&rot, a1, a2,  sina);
	mat_set(&rot, a2, a1, -sina);

	*angmom = mat_vec(&rot, angmom);

	/* transpose */
	mat_set(&rot, a1, a2, -sina);
	mat_set(&rot, a2, a1,  sina);

	*rotmat = mat_mat(rotmat, &rot);
}

/*
 * Rotation algorithm reference:
 *
 * Andreas Dullweber, Benedict Leimkuhler, and Robert McLachlan
 *
 * Symplectic splitting methods for rigid body molecular dynamics
 *
 * J. Chem. Phys. 107, 5840 (1997)
 */
static void rotate_body(struct body *body, double dt)
{
	double angle;

	/* rotate about x axis */
	angle = 0.5 * dt * body->angmom.x * body->inertia_inv.x;
	rotate_step(1, 2, angle, &body->angmom, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y * body->inertia_inv.y;
	rotate_step(2, 0, angle, &body->angmom, &body->rotmat);

	/* rotate about z axis */
	angle = dt * body->angmom.z * body->inertia_inv.z;
	rotate_step(0, 1, angle, &body->angmom, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y * body->inertia_inv.y;
	rotate_step(2, 0, angle, &body->angmom, &body->rotmat);

	/* rotate about x axis */
	angle = 0.5 * dt * body->angmom.x * body->inertia_inv.x;
	rotate_step(1, 2, angle, &body->angmom, &body->rotmat);
}

static void update_step_nve(struct md *md)
{
	double dt = cfg_get_double(md->state->cfg, "time_step");

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	compute_forces(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;
	}
}

/*
 * NVT with Nose-Hoover thermostat:
 *
 * William G. Hoover
 *
 * Canonical dynamics: Equilibrium phase-space distributions
 *
 * Phys. Rev. A 31, 1695 (1985)
 */
static void update_step_nvt(struct md *md)
{
	struct nvt_data *data = (struct nvt_data *)md->data;

	double dt = cfg_get_double(md->state->cfg, "time_step");
	double target = cfg_get_double(md->state->cfg, "temperature");
	double tau = cfg_get_double(md->state->cfg, "thermostat_tau");

	double t0 = get_temperature(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * data->chi);
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * data->chi);
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * data->chi);

		body->angmom.x += 0.5 * dt * (body->torque.x -
					body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
					body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
					body->angmom.z * data->chi);

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	data->chi += 0.5 * dt * (t0 / target - 1.0) / tau / tau;
	data->chi_dt += 0.5 * dt * data->chi;

	compute_forces(md);

	double chi_init = data->chi;
	vec_t angmom_init[md->n_bodies], vel_init[md->n_bodies];

	for (size_t i = 0; i < md->n_bodies; i++) {
		angmom_init[i] = md->bodies[i].angmom;
		vel_init[i] = md->bodies[i].vel;
	}

	for (size_t iter = 1; iter <= MAX_ITER; iter++) {
		double chi_prev = data->chi;
		double ratio = get_temperature(md) / target;

		data->chi = chi_init + 0.5 * dt * (ratio - 1.0) / tau / tau;

		for (size_t i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;

			body->vel.x = vel_init[i].x + 0.5 * dt *
				(body->force.x / body->mass - vel_init[i].x * data->chi);
			body->vel.y = vel_init[i].y + 0.5 * dt *
				(body->force.y / body->mass - vel_init[i].y * data->chi);
			body->vel.z = vel_init[i].z + 0.5 * dt *
				(body->force.z / body->mass - vel_init[i].z * data->chi);

			body->angmom.x = angmom_init[i].x + 0.5 * dt *
				(body->torque.x - angmom_init[i].x * data->chi);
			body->angmom.y = angmom_init[i].y + 0.5 * dt *
				(body->torque.y - angmom_init[i].y * data->chi);
			body->angmom.z = angmom_init[i].z + 0.5 * dt *
				(body->torque.z - angmom_init[i].z * data->chi);
		}

		if (fabs(data->chi - chi_prev) < EPSILON)
			break;

		if (iter == MAX_ITER)
			msg("WARNING: NVT UPDATE DID NOT CONVERGE\n\n");
	}

	data->chi_dt += 0.5 * dt * data->chi;
}

/*
 * Reference
 *
 * Simone Melchionna, Giovanni Ciccotti, Brad Lee Holian
 *
 * Hoover NPT dynamics for systems varying in shape and size
 *
 * Mol. Phys. 78, 533 (1993)
 */
static void update_step_npt(struct md *md)
{
	struct npt_data *data = (struct npt_data *)md->data;

	double dt = cfg_get_double(md->state->cfg, "time_step");
	double t_tau = cfg_get_double(md->state->cfg, "thermostat_tau");
	double t_target = cfg_get_double(md->state->cfg, "temperature");
	double p_tau = cfg_get_double(md->state->cfg, "barostat_tau");
	double p_target = cfg_get_double(md->state->cfg, "pressure");

	double t_tau2 = t_tau * t_tau;
	double p_tau2 = p_tau * p_tau;
	double kbt = BOLTZMANN * t_target;

	double t0 = get_temperature(md);
	double p0 = get_pressure(md);
	double v0 = get_volume(md);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * (data->chi + data->eta));
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * (data->chi + data->eta));
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * (data->chi + data->eta));

		body->angmom.x += 0.5 * dt * (body->torque.x -
						body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
						body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
						body->angmom.z * data->chi);

		rotate_body(body, dt);
	}

	data->chi += 0.5 * dt * (t0 / t_target - 1.0) / t_tau2;
	data->chi_dt += 0.5 * dt * data->chi;
	data->eta += 0.5 * dt * v0 * (p0 - p_target) / md->n_bodies / kbt / p_tau2;

	vec_t com = get_system_com(md);
	vec_t pos_init[md->n_bodies];

	for (size_t i = 0; i < md->n_bodies; i++)
		pos_init[i] = md->bodies[i].pos;

	for (size_t iter = 1; iter <= MAX_ITER; iter++) {
		bool done = true;

		for (size_t i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;
			vec_t pos = wrap(md, &body->pos);

			vec_t v = {
				data->eta * (pos.x - com.x),
				data->eta * (pos.y - com.y),
				data->eta * (pos.z - com.z)
			};

			vec_t new_pos = {
				pos_init[i].x + dt * (body->vel.x + v.x),
				pos_init[i].y + dt * (body->vel.y + v.y),
				pos_init[i].z + dt * (body->vel.z + v.z)
			};

			done = done && vec_dist(&body->pos, &new_pos) < EPSILON;
			body->pos = new_pos;
		}

		if (done)
			break;

		if (iter == MAX_ITER)
			msg("WARNING: NPT UPDATE DID NOT CONVERGE\n\n");
	}

	vec_scale(&md->box, exp(dt * data->eta));
	check_fail(efp_set_periodic_box(md->state->efp, md->box.x, md->box.y, md->box.z));

	compute_forces(md);

	double chi_init = data->chi, eta_init = data->eta;
	vec_t angmom_init[md->n_bodies], vel_init[md->n_bodies];

	for (size_t i = 0; i < md->n_bodies; i++) {
		angmom_init[i] = md->bodies[i].angmom;
		vel_init[i] = md->bodies[i].vel;
	}

	for (size_t iter = 1; iter <= MAX_ITER; iter++) {
		double chi_prev = data->chi;
		double eta_prev = data->eta;
		double t_cur = get_temperature(md);
		double p_cur = get_pressure(md);
		double v_cur = get_volume(md);

		data->chi = chi_init + 0.5 * dt * (t_cur / t_target - 1.0) / t_tau2;
		data->eta = eta_init + 0.5 * dt * v_cur * (p_cur - p_target) /
							md->n_bodies / kbt / p_tau2;

		for (size_t i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;

			body->vel.x = vel_init[i].x + 0.5 * dt *
				(body->force.x / body->mass -
					vel_init[i].x * (data->chi + data->eta));
			body->vel.y = vel_init[i].y + 0.5 * dt *
				(body->force.y / body->mass -
					vel_init[i].y * (data->chi + data->eta));
			body->vel.z = vel_init[i].z + 0.5 * dt *
				(body->force.z / body->mass -
					vel_init[i].z * (data->chi + data->eta));

			body->angmom.x = angmom_init[i].x + 0.5 * dt *
				(body->torque.x - angmom_init[i].x * data->chi);
			body->angmom.y = angmom_init[i].y + 0.5 * dt *
				(body->torque.y - angmom_init[i].y * data->chi);
			body->angmom.z = angmom_init[i].z + 0.5 * dt *
				(body->torque.z - angmom_init[i].z * data->chi);
		}

		if (fabs(data->chi - chi_prev) < EPSILON &&
		    fabs(data->eta - eta_prev) < EPSILON)
			break;

		if (iter == MAX_ITER)
			msg("WARNING: NPT UPDATE DID NOT CONVERGE\n\n");
	}

	data->chi_dt += 0.5 * dt * data->chi;
}

static void print_info(const struct md *md)
{
	double energy = md->potential_energy;
	double invariant = md->get_invariant(md);
	double temperature = get_temperature(md);

	msg("%30s %16.10lf\n", "POTENTIAL ENERGY", energy);
	msg("%30s %16.10lf\n", "INVARIANT", invariant);
	msg("%30s %16.10lf\n", "TEMPERATURE", temperature);

	if (cfg_get_enum(md->state->cfg, "ensemble") == ENSEMBLE_TYPE_NPT) {
		double pressure = get_pressure(md) / BAR_TO_AU;

		msg("%30s %16.10lf\n", "PRESSURE", pressure);
	}

	if (cfg_get_bool(md->state->cfg, "enable_pbc")) {
		double x = md->box.x * BOHR_RADIUS;
		double y = md->box.y * BOHR_RADIUS;
		double z = md->box.z * BOHR_RADIUS;

		msg("%30s %9.3lf %9.3lf %9.3lf\n", "PERIODIC BOX SIZE", x, y, z);
	}

	msg("\n\n");
}

static void print_restart(const struct md *md)
{
	msg("    RESTART DATA\n\n");

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		char name[64];
		check_fail(efp_get_frag_name(md->state->efp, i, sizeof(name), name));

		double xyzabc[6] = { body->pos.x * BOHR_RADIUS,
				     body->pos.y * BOHR_RADIUS,
				     body->pos.z * BOHR_RADIUS };

		matrix_to_euler(&body->rotmat, xyzabc + 3, xyzabc + 4, xyzabc + 5);

		double vel[6] = { body->vel.x,
				  body->vel.y,
				  body->vel.z,
				  body->angmom.x * body->inertia_inv.x,
				  body->angmom.y * body->inertia_inv.y,
				  body->angmom.z * body->inertia_inv.z };

		print_fragment(name, xyzabc, vel);
	}

	msg("\n");
}

static struct md *md_create(struct state *state)
{
	struct md *md = xcalloc(1, sizeof(struct md));

	md->state = state;
	md->box = box_from_str(cfg_get_string(state->cfg, "periodic_box"));

	switch (cfg_get_enum(state->cfg, "ensemble")) {
		case ENSEMBLE_TYPE_NVE:
			md->get_invariant = get_invariant_nve;
			md->update_step = update_step_nve;
			break;
		case ENSEMBLE_TYPE_NVT:
			md->get_invariant = get_invariant_nvt;
			md->update_step = update_step_nvt;
			md->data = xcalloc(1, sizeof(struct nvt_data));
			break;
		case ENSEMBLE_TYPE_NPT:
			md->get_invariant = get_invariant_npt;
			md->update_step = update_step_npt;
			md->data = xcalloc(1, sizeof(struct npt_data));
			break;
		default:
			assert(0);
	}

	md->n_bodies = state->sys->n_frags;
	md->bodies = xcalloc(md->n_bodies, sizeof(struct body));

	double coord[6 * md->n_bodies];
	check_fail(efp_get_coordinates(state->efp, coord));

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->pos.x = coord[6 * i + 0];
		body->pos.y = coord[6 * i + 1];
		body->pos.z = coord[6 * i + 2];

		double a = coord[6 * i + 3];
		double b = coord[6 * i + 4];
		double c = coord[6 * i + 5];

		euler_to_matrix(a, b, c, &body->rotmat);

		body->vel.x = md->state->sys->frags[i].vel[0];
		body->vel.y = md->state->sys->frags[i].vel[1];
		body->vel.z = md->state->sys->frags[i].vel[2];

		set_body_mass_and_inertia(state->efp, i, body);

		body->angmom.x = md->state->sys->frags[i].vel[3] * body->inertia.x;
		body->angmom.y = md->state->sys->frags[i].vel[4] * body->inertia.y;
		body->angmom.z = md->state->sys->frags[i].vel[5] * body->inertia.z;

		md->n_freedom += 3;

		if (body->inertia.x > EPSILON)
			md->n_freedom++;
		if (body->inertia.y > EPSILON)
			md->n_freedom++;
		if (body->inertia.z > EPSILON)
			md->n_freedom++;
	}

	return (md);
}

static void velocitize(struct md *md)
{
	rand_init();

	double temperature = cfg_get_double(md->state->cfg, "temperature");
	double ke = temperature * BOLTZMANN * md->n_freedom / (2.0 * 6.0 * md->n_bodies);

	for (size_t i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		double vel = sqrt(2.0 * ke / body->mass);

		body->vel.x = vel * rand_normal();
		body->vel.y = vel * rand_normal();
		body->vel.z = vel * rand_normal();

		body->angmom.x = sqrt(2.0 * ke * body->inertia.x) * rand_normal();
		body->angmom.y = sqrt(2.0 * ke * body->inertia.y) * rand_normal();
		body->angmom.z = sqrt(2.0 * ke * body->inertia.z) * rand_normal();
	}
}

static void print_status(const struct md *md)
{
	print_geometry(md->state->efp);
	print_restart(md);
	print_info(md);

	fflush(stdout);
}

static void md_shutdown(struct md *md)
{
	free(md->bodies);
	free(md->data);
	free(md);
}

void sim_md(struct state *state)
{
	msg("MOLECULAR DYNAMICS JOB\n\n\n");

	struct md *md = md_create(state);

	if (cfg_get_bool(state->cfg, "velocitize"))
		velocitize(md);

	remove_system_drift(md);
	compute_forces(md);

	msg("    INITIAL STATE\n\n");
	print_status(md);

	for (int i = 1; i <= cfg_get_int(state->cfg, "max_steps"); i++) {
		md->update_step(md);

		if (i % cfg_get_int(state->cfg, "print_step") == 0) {
			msg("    STATE AFTER %d STEPS\n\n", i);
			print_status(md);
		}
	}

	md_shutdown(md);

	msg("MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n");
}
