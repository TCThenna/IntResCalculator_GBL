   6   -1
   0   0.05  93.66  1  0.00342 
 150   0.05  93.66  1  0.00342 
 300   0.05  93.66  1  0.00342 
 450   0.05  93.66  1  0.00342 
 600   0.05  93.66  1  0.00342 
 750   0.05  93.66  1  0.00342 

======================================================
Description of geometry input file for AnaTel

First line:  
  _nTelPlanes _iDUT 
          _nTelPlanes - number of active and passive planes in the setup
          _iDUT       - position of the DUT (0..._nTelPlanes-1), -1 for no DUT

Following _nTelPlanes lines:
  _planePosition _planeThickness _planeX0 _isActive _planeResolution
where:   _planePosition - plane position along the beam [mm]
         _planeThickness - plane thickness [mm]
         _planeX0        - radiation lenght in plane material [mm]
         _isActive       - flag for active planes (0/1)
         _planeResolution - active plane resolution

 
