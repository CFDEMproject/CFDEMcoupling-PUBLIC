
// set force and mass source at cell subCellID
{
    // locate subPoint
    vector subPosition = position + offset;
    label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
               
    //NP fprintf(lmp->screen,"cellID=%d, partCellId=%d\n",static_cast<int>(cellID),static_cast<int>(partCellId));

    if(partCellId!=cellID)
    {
        nUnEqual++;
    }
    if(partCellId<0)
    {
        nNotFound++;
    }
    nTotal++;

    if (partCellId >= 0)  // subPoint is in domain
    {
        // update voidfraction for each particle read
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];
        scalar particleVolume = volume/static_cast<scalar>(nPoints);
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
        
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }

        cellsSet++; // inc counter of cells affected

        //====================================================//
        // add sub particle representation
        bool createNew = true;
        label storeInIndex=0;
        for(int i=0; i < cellsPerParticle_[index][0] ; i++)
        {
            if(partCellId == particleCloud_.cellIDs()[index][i]) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }
        }

        if(createNew)
        {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0]-1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;

        }

        particleWeights[index][storeInIndex] += 1/static_cast<scalar>(nPoints);
        particleVolumes[index][storeInIndex] += particleVolume;
        //====================================================//
    }
}
