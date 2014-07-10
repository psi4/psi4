#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga++.h"

GA::PGroup* GA::PGroup::pgMirror  = NULL;
GA::PGroup* GA::PGroup::pgDefault = NULL;
GA::PGroup* GA::PGroup::pgWorld   = NULL;

/**
 * Private -- never called.
 */
GA::PGroup::PGroup(void)
{
}

/**
 * Constructors and Destructor of PGroup
 */
GA::PGroup::PGroup(int *plist, int size) 
{
    mPHandle = GA_Pgroup_create(plist, size);
}

GA::PGroup::~PGroup() 
{
    GA_Pgroup_destroy(mPHandle);
}

/**
 * Pgroup Methods
 */
GA::PGroup*
GA::PGroup::getDefault() 
{
    if(pgDefault == NULL) 
    {
       pgDefault = new PGroup();
       pgDefault->mPHandle = GA_Pgroup_get_default();
    }
    return pgDefault;
}


GA::PGroup*
GA::PGroup::getMirror() 
{
    if(pgMirror == NULL) 
    {
       pgMirror = new PGroup();
       pgMirror->mPHandle = GA_Pgroup_get_mirror();
    }
    return pgMirror;    
}


GA::PGroup*
GA::PGroup::getWorld() 
{
    if(pgWorld == NULL) 
    {
       pgWorld = new PGroup();
       pgWorld->mPHandle = GA_Pgroup_get_world();
    }
    return pgWorld;
}


void
GA::PGroup::setDefault(GA::PGroup* p_handle) 
{
    pgDefault->mPHandle = p_handle->mPHandle;
}

void
GA::PGroup::brdcst(void* buf, int lenbuf, int root)
{
    GA_Pgroup_brdcst(mPHandle, buf, lenbuf, root);    
}

void
GA::PGroup::gop(double *buf, int n, char* op)
{
    GA_Pgroup_dgop(mPHandle, buf, n, op);                               
}

void
GA::PGroup::gop(long *buf, int n, char* op)
{
    GA_Pgroup_lgop(mPHandle, buf, n, op);
}

void
GA::PGroup::gop(int *buf, int n, char* op)
{
    GA_Pgroup_igop(mPHandle, buf, n, op);
}

void
GA::PGroup::gop(float *buf, int n, char* op)
{
    GA_Pgroup_fgop(mPHandle, buf, n, op);
}

int
GA::PGroup::nodeid() 
{
    return GA_Pgroup_nodeid(mPHandle);
}

int
GA::PGroup::nodes() 
{
    return GA_Pgroup_nnodes(mPHandle);
}

void
GA::PGroup::sync()
{
    GA_Pgroup_sync(mPHandle);
}

