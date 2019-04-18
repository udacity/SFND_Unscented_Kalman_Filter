#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "render/render.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

struct lmarker
{
	double x, y;
	lmarker(double setX, double setY)
		: x(setX), y(setY)
	{}

};

struct rmarker
{
	double rho, phi, rho_dot;
	rmarker(double setRho, double setPhi, double setRhoDot)
		: rho(setRho), phi(setPhi), rho_dot(setRhoDot)
	{}

};

class Tools {
	public:
	/**
	* Constructor.
	*/
	Tools();
	
	/**
	* Destructor.
	*/
	virtual ~Tools();
	
	// Members
	std::vector<VectorXd> estimations;
	std::vector<VectorXd> ground_truth;
	
	double noise(double stddev);
	lmarker lidarSense(Car& car, pcl::visualization::PCLVisualizer::Ptr& viewer, long long timestamp, bool visualize);
	rmarker radarSense(Car& car, Car ego, pcl::visualization::PCLVisualizer::Ptr& viewer, long long timestamp, bool visualize);
	void ukfResults(Car car, pcl::visualization::PCLVisualizer::Ptr& viewer, double time, int steps);
	/**
	* A helper method to calculate RMSE.
	*/
	VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);
	
};

#endif /* TOOLS_H_ */
