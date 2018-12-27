#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>



typedef std::array<double, 3> Vector;


Vector add(const Vector& target, const Vector& orig) {

  Vector ans = {0};

  ans[0] = target[0] + orig[0];
  ans[1] = target[1] + orig[1];
  ans[2] = target[2] + orig[2];

  return ans;
}

Vector sub(const Vector& target, const Vector& orig) {

  Vector ans = {0};

  ans[0] = target[0] - orig[0];
  ans[1] = target[1] - orig[1];
  ans[2] = target[2] - orig[2];

  return ans;
}

Vector sqr(const Vector& target) {

  Vector ans = {0};

  ans[0] = target[0] * target[0];
  ans[1] = target[1] * target[1];
  ans[2] = target[2] * target[2];

  return ans;
}

Vector scalar_mult(const Vector& target, double scalar) {

  Vector ans = target;

  ans[0] *= scalar;
  ans[1] *= scalar;
  ans[2] *= scalar;

  return ans;
}

double sum_vec(const Vector& target) {
  return target[0] + target[1] + target[2];
}

void print_vec(const Vector& target) {
  std::cout << target[0] << " " << target[1] << " " << target[2];
}

class Body {

public:

  std::string name;
  double mass;
  Vector position;

  std::vector< Vector > trajectory;

  Vector velocity = {0};
  Vector acceleration = {0};

  Body(std::string name, unsigned long long mass,
        const Vector& position, const Vector& velocity) {

    this -> name = name;
    this -> mass = mass;

    this->velocity = velocity;
    this->position = position;

  };

};



class Integrator {


private:

  int n_timesteps;

  Vector partial_step(const Vector& f, const Vector& df, double scale) {

    return Vector{
      f[0] + df[0] * n_timesteps * scale,
      f[1] + df[1] * n_timesteps * scale,
      f[2] + df[2] * n_timesteps * scale
    };

  };

  int accelerate(int this_index) {

    double G = 6.67428e-11;

    Body* target_body = &bodies[this_index];

    target_body->acceleration = {0};

    Vector velocity_update = {0};
    Vector location_update = {0};

    for (int i = 0; i < bodies.size(); ++i) {


      if (this_index == i)
        continue;

      Body* other_body = &bodies[i];

      Vector k1 = {0};
      Vector k2 = {0};
      Vector k3 = {0};
      Vector k4 = {0};


      Vector dist = sub(other_body->position, target_body->position);
      Vector sqred = sqr(dist);

      double r = sqrt( sum_vec(sqred) );

      //std::cout << "Dist: " << r << "\n";
      //std::cout << other_body->mass << "\n";

      auto tmp = G * other_body->mass / (r*r*r) * 1e20;

      //std::cout << "Tmp: " << tmp << "\n";

      k1 = scalar_mult(dist, tmp);

      velocity_update = partial_step(target_body->velocity, k1, 0.5);
      location_update = partial_step(target_body->position, velocity_update, 0.5);
      k2 = scalar_mult(sub(other_body->position, location_update), tmp);

      velocity_update = partial_step(target_body->velocity, k2, 0.5);
      location_update = partial_step(target_body->position, velocity_update, 0.5);
      k3 = scalar_mult(sub(other_body->position, location_update), tmp);

      velocity_update = partial_step(target_body->velocity, k3, 1);
      location_update = partial_step(target_body->position, velocity_update, 1);
      k4 = scalar_mult(sub(other_body->position, location_update), tmp);

      Vector temp_vec = add(k1, k4);
      temp_vec = add(temp_vec, scalar_mult(k2, 2));
      temp_vec = add(temp_vec, scalar_mult(k3, 2));

      temp_vec = scalar_mult(temp_vec, 0.17);

//      std::cout << "Temp vec"; print_vec(temp_vec); std::cout << "\n";
      target_body->acceleration = add(target_body->acceleration, temp_vec);
    }

    return 0;
  };

public:

  std::vector<Body> bodies;

  Integrator(std::vector<Body>& bodies, int n_timesteps) {
    this->bodies = bodies;
    this->n_timesteps = n_timesteps;
  };

  int update_coordinates() {

    for (int i = 0; i < this->bodies.size(); ++i)
      accelerate(i);

    for (int i = 0; i < this->bodies.size(); ++i) {
      Body* this_body = &this->bodies[i];

      this_body->velocity = add( this_body->velocity,
                scalar_mult(this_body->acceleration, n_timesteps));
      this_body->position = add( this_body->position,
                scalar_mult(this_body->velocity, n_timesteps));

    }

    return 0;
  };

};



int main() {

  std::vector<Body> bodies;

  //We devide each mass by 1e20

  bodies.push_back( Body("Sun", 2e10, {0, 0, 0}, {0, 0, 0}) );
  bodies.push_back( Body("Mercury", 3.285e3, {0, 5e10, 0}, {47000, 0, 0}) );
  bodies.push_back( Body("Venus", 4.8e4, {0, 1.1e11, 0}, {35000, 0, 0}) );
  bodies.push_back( Body("Earth", 6e4, {0, 1.5e11, 0}, {30000, 0, 0}) );
  bodies.push_back( Body("Mars", 2.4e4, {0, 2.2e11, 0}, {24000, 0, 0}) );
  bodies.push_back( Body("Jupiter", 1e8, {0, 7.7e11, 0}, {13000, 0, 0}) ); // Check speed here
  bodies.push_back( Body("Saturn", 5.7e6, {0, 1.4e12, 0}, {9000, 0, 0}) );
  bodies.push_back( Body("Uranus", 8.7e5, {0, 2.8e12, 0}, {6835, 0, 0}) );
  bodies.push_back( Body("Neptune", 1e6, {0, 4.5e12, 0}, {5477, 0, 0}) );
  bodies.push_back( Body("Pluto", 1.3e2, {0, 7.3e12, 0}, {4748, 0, 0}) );






  Integrator integr(bodies, 1000);

/*
  for (int i = 0; i < 10; ++i) {
    integr.update_coordinates();
    std::cout << "Earth: " << integr.bodies[1].position[0] << " " << integr.bodies[1].position[1] << " " << integr.bodies[1].position[2] << "\n";
    std::cout << "Acceleration: " << integr.bodies[1].acceleration[0] << " " << integr.bodies[1].acceleration[1] << " " << integr.bodies[1].acceleration[2] << "\n";
    std::cout << "Velocity: " << integr.bodies[1].velocity[0] << " " << integr.bodies[1].velocity[1] << " " << integr.bodies[1].velocity[2] << "\n";

  }
*/

  std::ofstream myfile;
  myfile.open("two_bodies.csv");

  std::string out;

  for (int i = 0; i < 80000000; ++i) {
    integr.update_coordinates();

    for (auto b : integr.bodies) {
      out = "";
      out += b.name + ",";
      out += std::to_string(b.position[0]) + ",";
      out += std::to_string(b.position[1]) + ",";
      out += std::to_string(b.position[2]) + "\n";

      myfile << out;
    }
  }

  return 0;
}
