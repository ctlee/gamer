{% set fullname = fullname | replace(module ~ ".", "") %}
{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% if methods %}
      .. rubric:: Methods

      .. autosummary::
      {% for item in methods %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}

   {% if attributes %}
      .. rubric:: Attributes

      .. autosummary::
      {% for item in attributes %}
         {{ name }}.{{ item }}
      {% endfor %}
   {% endif %}
