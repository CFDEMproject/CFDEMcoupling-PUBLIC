
  //========================================================================//
  scalar countCell=0;              // number of cells touched by particles
  int points=0;                    // number of particles and sub-points
  scalar totalParticleWeights=0;   // total weight of all particles and sub-points
  vector totalForce_array(0,0,0);  // total force on particles based on particle array
  vector totalForce_field(0,0,0);  // forceField of forceM(), used to calc Ksl
  scalar particleVolume_radius=0;  // total particle voulme based on radius array
  scalar particleVolume_field=0;   // total particle voulme based on voidfractionfield
  scalar particleVolume_array=0;   // total particle voulme based on particle array
  scalar meanR_array=0;            // mean particle radius based on particle array
  vector meanUs_array(0,0,0);      // mean solid velocity based on particle array
  vector meanUs_field(0,0,0);      // mean solid velocity based on field Us
  scalar meanAlpha_field=0;        // mean voidfraction
  vector meanU_field(0,0,0);       // mean fluid velocity based on field U
  //========================================================================//

  // get particle based data
  for(int index = 0;index <  numberOfParticles(); ++index)
  {
      for(int i=0;i<3;i++){
          totalForce_array[i] += impForces_[index][i]; // in giveDEMdata() exp imp and DEM are summed
          meanUs_array[i] += velocities_[index][i];
      }
      meanR_array += radii_[index][0];
      particleVolume_radius += radii_[index][0]*radii_[index][0]*radii_[index][0]*4./3.*3.1415;
      
      // loop subCells
      for(int subCell=0;subCell<voidFractionM().cellsPerParticle()[index][0];subCell++)
      {
          points++;
          totalParticleWeights += particleWeights_[index][subCell];
          particleVolume_array += particleVolumes_[index][subCell];
      }
  }

  // get field based data
  forAll(alpha,cellI)
  {
      // check if a particle is inside cell
      bool particleInside=false;
      for(int index = 0;index <  numberOfParticles(); ++index){
          for(int subCell=0;subCell<voidFractionM().cellsPerParticle()[index][0];subCell++){
              if(cellIDs_[index][subCell] == cellI){
                  particleInside=true;
                  break;
              }
          }
      }

      if(particleInside)
      {
          countCell++;
          meanAlpha_field += alpha[cellI];
          meanU_field += U[cellI];
          meanUs_field += Us[cellI];
          particleVolume_field += (1-alpha[cellI])*alpha.mesh().V()[cellI];
          totalForce_field += forceM(0).expParticleForces()[cellI]+forceM(0).impParticleForces()[cellI];
      }
  }

  // averaging
  if(countCell>0)
  {
      meanAlpha_field /= countCell;
      meanU_field /= countCell; 
      meanUs_field /= countCell;
  }
  else
  {
      meanAlpha_field = 0;
      meanU_field = vector(0,0,0); 
      meanUs_field = vector(0,0,0);
  }
  meanUs_array /= numberOfParticles()+SMALL;   
  meanR_array /= numberOfParticles()+SMALL;   

  Info <<"=============================================================================" << endl;
  Info << "Debug Info, only serial and not tested!" << endl;
  Info <<"  numberOfParticles_ = "<< numberOfParticles() << " != " << endl;
  Info <<"totalParticleWeights = "<< totalParticleWeights << endl;
  Info <<"   points= "<< points << endl;
  Info <<"countCell= "<< countCell << endl;
  Info <<"                totalForce_array = "<< mag(totalForce_array)<< " != " << endl;
  Info <<"                totalForce_field = "<< mag(totalForce_field) << endl;
  Info <<"            particleVolume_field = "<< particleVolume_field << " != " << endl;
  Info <<"            particleVolume_array = "<< particleVolume_array << " != " << endl;
  Info <<"           particleVolume_radius = "<< particleVolume_radius << endl;
  Info <<"meanUs_field = "<< mag(meanUs_field) << " ~= " << endl;
  Info <<"meanUs_array = "<< mag(meanUs_array) << endl;
  Info <<"meanU_field = "<< mag(meanU_field) << endl;
  Info <<"meanAlpha_field = "<< meanAlpha_field << endl;
  Info <<"meanR_array = "<< meanR_array << endl;
  Info <<"=============================================================================" << endl;
  Info << endl;
