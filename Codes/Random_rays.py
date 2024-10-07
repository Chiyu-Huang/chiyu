import numpy as np

def to_cl(kpc):
    return kpc*3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22) #have to change accordingly with the specific file

def start_end(k_unnorm, b = None, maxb = 25, minb = 5):
    '''
    Takes a input direction of light ray, then outputs the starting and 
    ending coordinates of the light ray that spans exactly the entire box.

    Args:
        k_unnorm (list) : direction of the desired light ray (not normalised)

    Returns:
        list : Starting point of the ray
        list : Ending point of the ray
    
    Raises:
        Exception : 'k must have 3 components (x, y, z)'
    '''
    start = []
    end = []
    while True:
        if len(k_unnorm) != 3:
            raise Exception('k must have 3 components (x, y, z)')
        kx, ky, kz = k_unnorm
        k = k_unnorm/np.sqrt(kx*kx + ky*ky + kz*kz)

        #   Define the orthonormal axis on the plane perpendicular to k, cutting
        #   through the centre of the galaxy
        x_hat = [1, 0, 0] 
        x_plane = np.cross(k, x_hat)
        x_plane = x_plane/np.sqrt(x_plane[0]**2 + x_plane[1]**2 + x_plane[2]**2)
        y_hat = [0, 1, 0]
        y_plane = np.cross(k, y_hat)
        y_plane = y_plane/np.sqrt(y_plane[0]**2 + y_plane[1]**2 + y_plane[2]**2)

        #   Random point on that plane, used to define the ray as it will 
        #   pass this point. Then also find when it intersects the box.
        impact_param = to_cl(np.random.uniform(minb, maxb))
        # print(impact_param)
        phi = np.random.uniform(0, 2*np.pi)
        x = impact_param*np.cos(phi)
        y = impact_param*np.sin(phi)

        # Find where ther boundaries of the box is reached
        rand_point = x*x_plane + y*y_plane + np.array([0.5, 0.5, 0.5])
        rand_point = np.array(rand_point)

        t0 = np.where(k > 0, -rand_point / k, np.where(k < 0, (1 - rand_point)/k, -1e10))
        t1 = np.where(k > 0, (1-rand_point) / k, np.where(k < 0, -rand_point/k, 1e10))
        t0 = np.max(t0)
        t1 = np.min(t1)
        # print(t0, t1)



        start = t0*k + rand_point
        end = t1*k + rand_point
        th = 1e-14
        if np.min(start)<0 or np.max(start)>1 or np.min(end)<0 or np.max(end)>1:
            if is_close(np.min(start), 0, th):
                start[np.argmin(start)]=0
            if is_close(np.max(start), 1, th):
                start[np.argmax(start)] = 1
            if is_close(np.min(end), 0, th):
                end[np.argmin(end)] = 0
            if is_close(np.max(end), 1, th):
                end[np.argmax(end)] = 1
            if np.min(end)<1e-14 or np.min(start)<1e-14 or np.max(end)>1+1e-14 or np.max(start)>1+1e-14:
                # print(start, end)
                continue
        else:
            break
    # print(start, end)
    return (start, end, np.sqrt(x*x + y*y))

def is_close(a, b, th):
    return np.abs(a-b)<th