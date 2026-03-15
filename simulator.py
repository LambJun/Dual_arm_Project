import pybullet 
import pybullet_data
import numpy as np
import numpy.matlib
import math 
import pybullet_utils
import time

from pybullet_utils import bullet_client as bc
from pybullet_utils import urdfEditor as ed

urdf_path = 'C:\\Users\\David\\OneDrive - National University of Singapore\\PostDoc_NUS\\Code\\urdf\\sw2urdf_test\\urdf_file\\sw2urdf_test';
urdf_path2 = 'C:\\Users\\David\\OneDrive - National University of Singapore\\PhD_work\\Code\\TAMP_HighLevel\\Visual_Studio\\PythonGrammarTAMP';

# Simulate the robot, the environment, and the interactions
class Simulator:

    # Class constructor
    def __init__(self):
        self.object1 = [[]];
        self.object2 = [[]];
        self.available_objects = [];
        self.attached_objects = [];
        self.pitcher = 0;
        self.mug = 0;

    ##  Convert numpy values from millimetres to metres.
    #   @param in_metres Array of values to convert in millimitres.
    #   @return Numpy array of same dimensions as in_metres with all values in metres
    def toMilli(self,in_metres):
        # Check that the coords vector is a column vector
        if np.ndim(in_metres)!=2 or np.size(in_metres,1)!=1 or np.size(in_metres,0)!=3:
            raise Exception("[simulator.py] The vector of coordinates is supposed to be a column vector with three rows (x,y,z). Here it is not.");
        # Check that the length between the simulation and robot's world coordinates origins is not negative or equal to zero
        return np.multiply(in_metres,0.001);#element-wise multiplication

    ##	Load and execute the pybullet simulator with the set of joint positions
    #	@param in_degrees Numpy matrix of dimensions 1x5 containing the joint positions (in radiians) for the left arm's joints
    #   @return numpy matrix containing values in degrees.
    def toDegrees(self,in_radians):
        return np.multiply(in_radians,180/math.pi);

    ##	Load and execute the pybullet simulator with the set of joint positions
    #	@param in_degrees Numpy matrix of dimensions 1x5 containing the joint positions (in degrees) for the left arm's joints
    #   @return numpy matrix containing values in radians.
    def toRadians(self,in_degrees):
        return np.multiply(in_degrees,math.pi/180);

    ##	Transform the coordinates from robot's world coordinates into simulation's world coordinates
    #	@param l0 Normal distance between the simulation and robot's world coordinates origins. In millimetres.
    #   @param coords Coordinates of the point in the robot's world frame. Numpy array. In millimetres.
    #   @return numpy array containing three coordinates: (x,y,z) expressed in the simulation world frame
    def transform2Sim(self,l0,coords):
        # Check that the coords vector is a column vector
        if np.ndim(coords)!=2 or np.size(coords,1)!=1 or np.size(coords,0)!=3:
            raise Exception("[simulator.py] The vector of coordinates is supposed to be a column vector with three rows (x,y,z). Here it is not.");
        # Check that the length between the simulation and robot's world coordinates origins is not negative or equal to zero
        if l0<=0:
            raise Exception("[simulator.py] The length cannot be negative or equal to zero. Here it is.");
        # Transformation matrix between the simulation and robot's world frames and multiplication.
        T = np.array([[0,-1,0,0],[1,0,0,0],[0,0,1,l0],[0,0,0,1]]);
        # Append one to the end of the coordinates vector. Otherwise, the arrays multiplication will raise an error.
        coords = np.matmul(T,np.vstack((coords,np.array([1]))));
        return coords[0:3,:];

    ##	Load and execute the pybullet simulator with the set of joint positions. Be aware that the two arms will move to the back for the initialization
    #	@input target_joint_positions Python list containing the joint positions (in radians)
    #   @input threshold Float. The threshold down below the object is considered attached to the end-effector (in metres)
    def loadSimulator(self,target_joint_positions,object1,object2,threshold):

        object1 = self.toMilli(self.transform2Sim(1120,object1[0:3,:]));# Position of object1 is given with respect to the robot's frame. Then, it is converted into the World coordinates
        object2 = self.toMilli(self.transform2Sim(1120,object2[0:3,:]));# Position of object2 is given with respect to the robot's frame. Then, it is converted into the World coordinates
        physicsClient = pybullet.connect(pybullet.GUI);#or pybullet.DIRECT for non-graphical version
        pybullet.setAdditionalSearchPath(pybullet_data.getDataPath()); #used by loadURDF
        #pybullet.setGravity(0,0,-10);
        pybullet.setRealTimeSimulation(1);
        pybullet.setTimeStep(0.000001);#In seconds
        planeId = pybullet.loadURDF("plane.urdf");
        birobot = pybullet.loadURDF(urdf_path+"\\urdf\\sw2urdf_test.urdf",[0,0,1.12],pybullet.getQuaternionFromEuler([0,0,0]));# First ensemble is position. Second ensemble is orientation.
        jointIndices = np.arange(1,pybullet.getNumJoints(birobot),1,dtype=int);

        ed0 = ed.UrdfEditor();
        ed0.initializeFromBulletBody(birobot,physicsClient);

        time.sleep(10);

        while 1:
            pybullet.stepSimulation();

        pybullet.disconnect();

    ##	Run the pybullet simulator
    #	@input angles_left Numpy matrix of dimensions 1x5 containing the joint positions (in degrees) for the left arm's joints
    #	@input angles_right Numpy matrix of dimensions 1x5 containing the joint positions (in degrees) for the right arm's joints
    #   @input threshold The threshold down below the object is considered to be attached to the end-effector
    def runSimulator(self,angles_left,angles_right,object1,object2,threshold):
        angles_left = self.toRadians(angles_left);
        angles_right = self.toRadians(angles_right);
        angles_left_sim = np.zeros((np.size(angles_left,0),np.size(angles_left,1)));
        angles_right_sim = np.zeros((np.size(angles_right,0),np.size(angles_right,1)));

        # The angles passed to runSimulator are all positive. However, the simulator's directions are different. Therefore, a mapping must take place.
        angles_left = np.multiply(angles_left, np.matlib.repmat(np.array([-1,-1,-1,-1,-1]),np.size(angles_left,0),1));
        angles_right = np.multiply(angles_right, np.matlib.repmat(np.array([1,-1,1,-1,1]),np.size(angles_right,0),1));
        angles_left_sim[:,0] = angles_left[:,0];
        angles_left_sim[:,1] = angles_left[:,1];
        angles_left_sim[:,2] = angles_left[:,2];
        angles_left_sim[:,3] = angles_left[:,3];
        angles_left_sim[:,4] = angles_left[:,4];
        angles_right_sim[:,0] = angles_right[:,0];
        angles_right_sim[:,1] = angles_right[:,1];
        angles_right_sim[:,2] = angles_right[:,2];
        angles_right_sim[:,3] = angles_right[:,3];
        angles_right_sim[:,4] = angles_right[:,4];
        target_joint_positions = list(np.hstack((angles_right_sim,angles_left_sim)));
        self.loadSimulator(target_joint_positions,object1,object2,threshold);

if __name__ == "__main__":
    sim = Simulator();
    angles_left = np.matlib.repmat(np.array([0,0,0,0,0]),100000,1);
    angles_right = np.matlib.repmat(np.array([0,0,0,0,0]),100000,1);
    object1 = [[181.0709334,-542.0914752,-260,0.02527894,0.01340268,-0.64181416,-0.76632625]];#pitcher at Pst: 5. Be careful the translation vector must be in mm!
    object2 = [[50.4252055,-545.6029044,-340,0.00091483,-0.04028703,0.63988639,0.76741223]];# cup at Pst: 6. Be careful the translation vector must be in mm!
    object1 = np.transpose(np.array(object1,ndmin=2));
    object2 = np.transpose(np.array(object2,ndmin=2));
    sim.runSimulator(angles_left,angles_right,object1,object2,0.1);