classdef MasterOrbit 
    
  properties 
      Value  {oe_vector};
  end
  methods 
      function target_orbit = DesiredOrbit(oe_vector)
          
          target_orbit = desired_orbit(oe_vector,delta_oe);
          
      end
  end
end


