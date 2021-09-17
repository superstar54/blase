import bpy
from bpy.props import StringProperty, IntProperty, FloatProperty, CollectionProperty, BoolProperty
from bpy.types import PropertyGroup, UIList, Operator, Panel

# https://blender.stackexchange.com/questions/47840/is-bpy-props-able-to-create-a-list-of-lists

class ListItem(PropertyGroup):
    """Group of properties representing an item in the list."""

    name: StringProperty(name="name", description="Species", default="")
    species1: StringProperty(name="species1", description="Species", default="")
    species2: StringProperty(name="species2", description="Species", default="")
    bondlength: FloatProperty(name="bondlength", description="")
    polyhedra: BoolProperty(name="polyhedra", description="")
    search: BoolProperty(name="search", description="")


class MY_UL_List(UIList):
    """Demo UIList."""

    # Filter by the value of bondlength
    filter_by_bondlength: StringProperty(default='')

    # Invert the random property filter
    invert_filter_by_random: BoolProperty(default=False)

    # Order by random prop
    order_by_bondlength: BoolProperty(default=False)


    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):

        # We could write some code to decide which icon to use here...
        custom_icon = 'OBJECT_DATAMODE'

        # Make sure your code supports all 3 layout types
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text='%s_%s'%(item.species1, item.species2), icon = "COLOR")
            layout.prop(item, "species1")
            layout.prop(item, "species2")
            layout.prop(item, "bondlength")
            layout.prop(item, "polyhedra")
            layout.prop(item, "polyhedra")

        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text='', icon = custom_icon)


    def draw_filter(self, context, layout):
        """UI code for the filtering/sorting/search area."""

        layout.separator()
        col = layout.column(align=True)

        row = col.row(align=True)
        row.prop(self, 'filter_by_bondlength', text='', icon='VIEWZOOM')
        row.prop(self, 'invert_filter_by_random', text='', icon='ARROW_LEFTRIGHT')


    def filter_items(self, context, data, propname):
        """Filter and order items in the list."""

        # We initialize filtered and ordered as empty lists. Notice that 
        # if all sorting and filtering is disabled, we will return
        # these empty. 

        filtered = []
        ordered = []
        items = getattr(data, propname)

        # Filter
        if self.filter_by_bondlength:

            # Initialize with all items visible
            filtered = [self.bitflag_filter_item] * len(items)

            for i, item in enumerate(items):
                if item.bondlength != self.filter_by_bondlength:
                    filtered[i] &= ~self.bitflag_filter_item


        # Invert the filter
        if filtered and self.invert_filter_by_random:
            show_flag = self.bitflag_filter_item & ~self.bitflag_filter_item

            for i, bitflag in enumerate(filtered):
                if bitflag == filter_flag:
                    filtered[i] = self.bitflag_filter_item
                else:
                    filtered[i] &= ~self.bitflag_filter_item


        # Order by the length of bondlength
        if self.order_by_bondlength:
            sort_items = bpy.types.UI_UL_list.helper_funcs.sort_items_helper
            ordered = sort_items(items, lambda i: len(i.bondlength), True)


        return filtered, ordered



class LIST_OT_NewItem(Operator):
    """Add a new item to the list."""

    bl_idname = "my_list.new_item"
    bl_label = "Add a new item"

    def execute(self, context):
        context.scene.my_list.add()

        return{'FINISHED'}


class LIST_OT_DeleteItem(Operator):
    """Delete the selected item from the list."""

    bl_idname = "my_list.delete_item"
    bl_label = "Deletes an item"

    @classmethod
    def poll(cls, context):
        return context.scene.my_list

    def execute(self, context):
        my_list = context.scene.my_list
        index = context.scene.list_index

        my_list.remove(index)
        context.scene.list_index = min(max(0, index - 1), len(my_list) - 1)

        return{'FINISHED'}


class LIST_OT_MoveItem(Operator):
    """Move an item in the list."""

    bl_idname = "my_list.move_item"
    bl_label = "Move an item in the list"

    direction: bpy.props.EnumProperty(items=(('UP', 'Up', ""),
                                              ('DOWN', 'Down', ""),))

    @classmethod
    def poll(cls, context):
        return context.scene.my_list

    def move_index(self):
        """ Move index of an item render queue while clamping it. """

        index = bpy.context.scene.list_index
        list_length = len(bpy.context.scene.my_list) - 1  # (index starts at 0)
        new_index = index + (-1 if self.direction == 'UP' else 1)

        bpy.context.scene.list_index = max(0, min(new_index, list_length))

    def execute(self, context):
        my_list = context.scene.my_list
        index = context.scene.list_index

        neighbor = index + (-1 if self.direction == 'UP' else 1)
        my_list.move(neighbor, index)
        self.move_index()

        return{'FINISHED'}


class PT_ListExample(Panel):
    """Demo panel for UI list Tutorial."""

    bl_label = "UI_List Demo"
    bl_idname = "SCENE_PT_LIST_DEMO"
    bl_space_type  = "VIEW_3D"
    # bl_space_type = 'PROPERTIES'
    # bl_region_type = 'WINDOW'
    bl_region_type = "UI"
    bl_category = "uilist"
    # bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        row = layout.row()
        row.template_list("MY_UL_List", "The_List", scene,
                          "my_list", scene, "list_index")

        row = layout.row()
        row.operator('my_list.new_item', text='NEW')
        row.operator('my_list.delete_item', text='REMOVE')
        row.operator('my_list.move_item', text='UP').direction = 'UP'
        row.operator('my_list.move_item', text='DOWN').direction = 'DOWN'

        if scene.list_index >= 0 and scene.my_list:
            item = scene.my_list[scene.list_index]

            layout.row().prop(item, 'name')
            layout.row().prop(item, 'species1')
            layout.row().prop(item, 'species2')
            layout.row().prop(item, 'bondlength')
            layout.row().prop(item, 'polyhedra')
            layout.row().prop(item, 'search')
