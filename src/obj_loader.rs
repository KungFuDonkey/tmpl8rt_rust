use crate::material::Material;
use crate::math::*;
use crate::objects::mesh::Mesh;

pub fn load_obj(path: &std::path::Path, obj_idx: usize, mat_idx: usize, transform: &Mat4) -> (Vec<Mesh>, Vec<Material>)
{
    let mut options = tobj::LoadOptions::default();
    options.triangulate = true;

    let (models, _) = tobj::load_obj(path, &options).expect("Failed to load obj file");

    let mut meshes: Vec<Mesh> = Vec::new();
    let mats: Vec<Material> = Vec::new();

    let mut obj_idx = obj_idx;
    //let mut mat_idx = mat_idx;

    for m in &models
    {
        let mesh = &m.mesh;
        let mat_id = mesh.material_id;
        /*let mut matx_id = 0;
        if mat_id != None
        {
            matx_id = mat_id.unwrap();
        }*/

        let mut vertices: Vec<Float3> = Vec::with_capacity(mesh.positions.len() / 3);
        let mut triangles: Vec<[usize;3]> = Vec::with_capacity(mesh.indices.len() / 3);

        for vtx in 0..mesh.positions.len() / 3
        {
            vertices.push(Float3::from_xyz(
                mesh.positions[3 * vtx + 0],
                mesh.positions[3 * vtx + 1],
                mesh.positions[3 * vtx + 2]
            ));
        }

        for vtx in 0..mesh.indices.len() / 3
        {
            triangles.push([
                mesh.indices[3 * vtx + 0] as usize,
                mesh.indices[3 * vtx + 1] as usize,
                mesh.indices[3 * vtx + 2] as usize,
            ])
        }

        meshes.push(Mesh::from_data(obj_idx, 0, *transform, vertices, triangles));
        obj_idx += 1;
    }


    return (meshes, mats);
}