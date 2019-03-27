module history_mod

  use params_mod
  use mesh_mod
  use io_mod
  use log_mod
  use string_mod
  use static_mod
  use state_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  interface history_write
    module procedure history_write_state
  end interface history_write

contains

  subroutine history_init()

    call io_init()

    if (output_file_prefix /= 'N/A') then
      call io_create_dataset(desc=case_name, file_prefix=output_file_prefix, frames_per_file=frames_per_file)
    else
      call io_create_dataset(desc=case_name, file_prefix=case_name, frames_per_file=frames_per_file)
    end if

    call io_add_meta('source',        'MCV_SW')
    call io_add_meta('dt',            dt)
    call io_add_meta('time_scheme',   time_scheme)
    call io_add_meta('author',        'N/A')
    call io_add_meta('on_a_sphere',   'YES')
    call io_add_meta('sphere_radius', radius)
    ! Dimensions
    call io_add_dim('Time',           add_var=.true.)
    call io_add_dim('nCells',         size=nCells)
    call io_add_dim('nEdges',         size=nEdges)
    call io_add_dim('nVertices',      size=nVertices)
    call io_add_dim('TWO',            size=2)
    call io_add_dim('vertexDegree',   size=vertexDegree)
    call io_add_dim('maxEdges',       size=maxEdges)
    call io_add_dim('maxEdges2',      size=maxEdges2)
    ! Mesh parameters
    call io_add_var('lonCell',        long_name='Longitude on the cell',                       units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('latCell',        long_name='Latitude on the cell',                        units='radian', dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('xCell',          long_name='Cartesian X on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('yCell',          long_name='Cartesian Y on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('zCell',          long_name='Cartesian Z on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('indexToCellID',  long_name='Global cell ID',                              units='1',      dim_names=['nCells      '],                 data_type='integer')
    call io_add_var('lonEdge',        long_name='Longitude on the edge',                       units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('latEdge',        long_name='Latitude on the edge',                        units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('xEdge',          long_name='Cartesian X on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('yEdge',          long_name='Cartesian Y on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('zEdge',          long_name='Cartesian Z on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('indexToEdgeID',  long_name='Global edge ID',                              units='1',      dim_names=['nEdges      '],                 data_type='integer')
    call io_add_var('lonVertex',      long_name='Longitude on the vertex',                     units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('latVertex',      long_name='Latitude on the vertex',                      units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('xVertex',        long_name='Cartesian X on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('yVertex',        long_name='Cartesian Y on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('zVertex',        long_name='Cartesian Z on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('indexToVertexID',long_name='Global vertex ID',                            units='1',      dim_names=['nVertices   '],                 data_type='integer')
    call io_add_var('areaCell',       long_name='Primary cell area',                           units='m2',     dim_names=['nCells      '],                 data_type='real(8)')
    call io_add_var('areaTriangle',   long_name='Dual cell area',                              units='m2',     dim_names=['nVertices   '],                 data_type='real(8)')
    call io_add_var('areaEdge',       long_name='Defined edge area',                           units='m2',     dim_names=['nEdges      '],                 data_type='real(8)')
    call io_add_var('nEdgesOnCell',   long_name='Edge number on the cell',                     units='1',      dim_names=['nCells      '],                 data_type='integer')
    call io_add_var('nEdgesOnEdge',   long_name='Edge number to reconstruct tangent velocity', units='1',      dim_names=['nEdges      '],                 data_type='integer')
    call io_add_var('cellsOnCell',    long_name='Cell indices that surround cell',             units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('cellsOnEdge',    long_name='Cell indices that saddle cell',               units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    call io_add_var('cellsOnVertex',  long_name='Cell indices that surround vertex',           units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call io_add_var('edgesOnCell',    long_name='Edge indices on the cell',                    units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('edgesOnEdge',    long_name='Edge indices to reconstruct tangent velocity',units='1',      dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer')
    call io_add_var('edgesOnVertex',  long_name='Edge indices on the vertex',                  units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
    call io_add_var('verticesOnCell', long_name='Vertex indices on the cell',                  units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
    call io_add_var('verticesOnEdge', long_name='Vertex indices on the edge',                  units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
    ! Dynamical variables
    call io_add_var('u',              long_name='Normal wind on the edge',                     units='m s-1',  dim_names=['nEdges   ', 'Time     '],       data_type='real(8)')
    call io_add_var('h',              long_name='Geopotential height on the cell',             units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('pv',             long_name='Potential vorticity on the cell',             units='',       dim_names=['nVertices', 'Time     '],       data_type='real(8)')
    call io_add_var('div',            long_name='Divergence',                                  units='s-1',    dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
    call io_add_var('tm',             long_name='total mass',                                  units='m2 s-2', dim_names=['Time     '],                    data_type='real(8)')
    call io_add_var('te',             long_name='total energy',                                units='m4 s-4', dim_names=['Time     '],                    data_type='real(8)')

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state(state, static)

    type(state_type),  intent(inout) :: state
    type(static_type), intent(in   ) :: static

    call div_operator(state%edge%u, state%cell%div)

    call io_start_output()
    call io_output('lonCell',         lonCell)
    call io_output('latCell',         latCell)
    call io_output('xCell',           xCell)
    call io_output('yCell',           yCell)
    call io_output('zCell',           zCell)
    call io_output('indexToCellID',   indexToCellID)
    call io_output('lonEdge',         lonEdge)
    call io_output('latEdge',         latEdge)
    call io_output('xEdge',           xEdge)
    call io_output('yEdge',           yEdge)
    call io_output('zEdge',           zEdge)
    call io_output('indexToEdgeID',   indexToEdgeID)
    call io_output('lonVertex',       lonVertex)
    call io_output('latVertex',       latVertex)
    call io_output('xVertex',         xVertex)
    call io_output('yVertex',         yVertex)
    call io_output('zVertex',         zVertex)
    call io_output('indexToVertexID', indexToVertexID)
    call io_output('areaCell',        areaCell)
    call io_output('areaTriangle',    areaTriangle)
    call io_output('areaEdge',        areaEdge)
    call io_output('nEdgesOnCell',    nEdgesOnCell)
    call io_output('nEdgesOnEdge',    nEdgesOnEdge)
    call io_output('cellsOnCell',     cellsOnCell)
    call io_output('cellsOnEdge',     cellsOnEdge)
    call io_output('cellsOnVertex',   cellsOnVertex)
    call io_output('edgesOnCell',     edgesOnCell)
    call io_output('edgesOnEdge',     edgesOnEdge)
    call io_output('edgesOnVertex',   edgesOnVertex)
    call io_output('verticesOnCell',  verticesOnCell)
    call io_output('verticesOnEdge',  verticesOnEdge)
    call io_output('u',               state%edge%u)
    call io_output('h',               (state%cell%gd + static%cell%ghs) / g)
    call io_output('pv',              state%vertex%pv)
    call io_output('div',             state%cell%div)
    call io_output('tm',              state%total_mass)
    call io_output('te',              state%total_energy)
    call io_end_output()

  end subroutine history_write_state

end module history_mod
