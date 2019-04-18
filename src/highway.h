/* \author Aaron Brown */
// Handle logic for creating traffic on highway and animating it

#include "render/render.h"
#include "tools.h"

class Highway
{
public:

	std::vector<Car> traffic;
	Car egoCar;
	Tools tools;
	std::vector<bool> trackCars;

	Highway(pcl::visualization::PCLVisualizer::Ptr& viewer)
	{

		// Set which cars to track with UKF
		trackCars = {true,true,true};

		tools = Tools();
	
		egoCar = Car(Vect3(0, 0, 0), Vect3(4, 2, 2), Color(0, 1, 0), 0, 0, 2, "egoCar");
		
		Car car1(Vect3(-10, 4, 0), Vect3(4, 2, 2), Color(0, 0, 1), 5, 0, 2, "car1");
		
		std::vector<accuation> car1_instructions;
		accuation a = accuation(0.5*1e6, 0.5, 0.0);
		car1_instructions.push_back(a);
		a = accuation(2.2*1e6, 0.0, -0.2);
		car1_instructions.push_back(a);
		a = accuation(3.3*1e6, 0.0, 0.2);
		car1_instructions.push_back(a);
		a = accuation(4.4*1e6, -2.0, 0.0);
		car1_instructions.push_back(a);
	
		car1.setInstructions(car1_instructions);
		if( trackCars[0] )
		{
			UKF ukf1;
			car1.setUKF(ukf1);
		}
		traffic.push_back(car1);
		
		Car car2(Vect3(25, -4, 0), Vect3(4, 2, 2), Color(0, 0, 1), -6, 0, 2, "car2");
		std::vector<accuation> car2_instructions;
		a = accuation(4.0*1e6, 3.0, 0.0);
		car2_instructions.push_back(a);
		a = accuation(8.0*1e6, 0.0, 0.0);
		car2_instructions.push_back(a);
		car2.setInstructions(car2_instructions);
		if( trackCars[1] )
		{
			UKF ukf2;
			car2.setUKF(ukf2);
		}
		traffic.push_back(car2);
	
		Car car3(Vect3(-12, 0, 0), Vect3(4, 2, 2), Color(0, 0, 1), 1, 0, 2, "car3");
		std::vector<accuation> car3_instructions;
		a = accuation(0.5*1e6, 2.0, 1.0);
		car3_instructions.push_back(a);
		a = accuation(1.0*1e6, 2.5, 0.0);
		car3_instructions.push_back(a);
		a = accuation(3.2*1e6, 0.0, -1.0);
		car3_instructions.push_back(a);
		a = accuation(3.3*1e6, 2.0, 0.0);
		car3_instructions.push_back(a);
		a = accuation(4.5*1e6, 0.0, 0.0);
		car3_instructions.push_back(a);
		a = accuation(5.5*1e6, -2.0, 0.0);
		car3_instructions.push_back(a);
		a = accuation(7.5*1e6, 0.0, 0.0);
		car3_instructions.push_back(a);
		car3.setInstructions(car3_instructions);
		if( trackCars[2] )
		{
			UKF ukf3;
			car3.setUKF(ukf3);
		}
		traffic.push_back(car3);
	
		// render environment
		renderHighway(0,viewer);
		egoCar.render(viewer);
		car1.render(viewer);
		car2.render(viewer);
		car3.render(viewer);
	}
	
	void stepHighway(double egoVelocity, long long timestamp, int frame_per_sec, pcl::visualization::PCLVisualizer::Ptr& viewer)
	{
	
		// render highway environment with poles
		renderHighway(egoVelocity*timestamp/1e6, viewer);
		egoCar.render(viewer);
		
		for (int i = 0; i < traffic.size(); i++)
		{
			traffic[i].move((double)1/frame_per_sec, timestamp);
			traffic[i].render(viewer);
			// Sense surrounding cars with lidar and radar
			if(trackCars[i])
			{
				VectorXd gt(4);
				gt << traffic[i].position.x, traffic[i].position.y, traffic[i].velocity*cos(traffic[i].angle), traffic[i].velocity*sin(traffic[i].angle);
				tools.ground_truth.push_back(gt);
				tools.lidarSense(traffic[i], viewer, timestamp);
				tools.radarSense(traffic[i], egoCar, viewer, timestamp);
				tools.ukfResults(traffic[i],viewer);
				VectorXd estimate(4);
				double v  = traffic[i].ukf.x_(2);
    			double yaw = traffic[i].ukf.x_(3);
    			double v1 = cos(yaw)*v;
    			double v2 = sin(yaw)*v;
				estimate << traffic[i].ukf.x_[0], traffic[i].ukf.x_[1], v1, v2;
				tools.estimations.push_back(estimate);
	
			}
		}
		viewer->addText("Accuracy - RMSE:", 30, 300, 20, 1, 1, 1, "rmse");
		VectorXd rmse = tools.CalculateRMSE(tools.estimations, tools.ground_truth);
		viewer->addText(std::to_string(rmse[0]), 30, 275, 20, 1, 1, 1, "rmse_x");
		viewer->addText(std::to_string(rmse[1]), 30, 250, 20, 1, 1, 1, "rmse_y");
		viewer->addText(std::to_string(rmse[2]), 30, 225, 20, 1, 1, 1, "rmse_vx");
		viewer->addText(std::to_string(rmse[3]), 30, 200, 20, 1, 1, 1, "rmse_vy");
		
	}
	
};